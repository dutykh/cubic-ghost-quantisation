# Authors: Davide Batic and Denys Dutykh (Mathematics Department, Khalifa University of
# Science and Technology, Abu Dhabi, UAE).
# Copyright (C) 2026 Davide Batic and Denys Dutykh.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this repository.
"""Finite compressions of the cubic ghost generator.

The generator is m (N + 1/2) + c B**3 on the occupation-number basis e_0, e_1, ...,
with B e_n = sqrt(n+1) e_{n+1} - sqrt(n) e_{n-1} and c = lambda / (6 (2m)**(3/2)).
The dimension of a compression is N=J+1: states e_0,...,e_{N-1}.  The cubic matrix
is P_J B**3 P_J, with the power formed before compression, and not the cube of
the compressed matrix.

The module holds the matrix in double precision, with the condition number and
the residual of each eigenvalue; the same matrix in Decimal arithmetic, with the
determinant and the secant iteration that give its real levels to many digits;
and a fourth-order finite-difference calculation in coordinate space for the
lowest levels of the imaginary cubic oscillator, which does not use the matrix.
All quantities are dimensionless when m=1 and lambda=3.
"""
from decimal import Decimal, localcontext

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def ghost_matrix(n, mass=1.0, coupling=3.0):
    """The exact finite compression, in double precision."""
    k = np.arange(1, n, dtype=float)
    cubic = np.diag(-3 * k**1.5, -1) + np.diag(3 * k**1.5, 1)
    if n >= 4:
        k = np.arange(1, n - 2, dtype=float)
        third = np.sqrt(k * (k + 1) * (k + 2))
        cubic += np.diag(third, -3) - np.diag(third, 3)
    coefficient = coupling / (6 * (2 * mass)**1.5)
    return mass * np.diag(np.arange(n) + 0.5) + coefficient * cubic


def spectrum_and_condition_numbers(n):
    """Eigenvalues of the compression of size n, sorted by real part, with the
    condition number and the relative right residual of each."""
    matrix = ghost_matrix(n)
    values, left, right = la.eig(matrix, left=True, right=True)
    order = np.lexsort((values.imag, values.real))
    values, left, right = values[order], left[:, order], right[:, order]
    overlaps = np.sum(left.conj() * right, axis=0)
    # scipy.linalg.eig returns left and right vectors of Euclidean norm one.
    condition = 1 / np.abs(overlaps)
    residual = np.linalg.norm(matrix @ right - right * values, axis=0)
    residual /= (la.norm(matrix, 2) + np.abs(values))
    return values, condition, residual


def decimal_matrix(n):
    """Seven-diagonal compression at m=1, lambda=3.

    Call inside the required decimal.localcontext.  Entries are evaluated
    in Decimal arithmetic from integers, including the square roots.
    """
    D = Decimal
    coefficient = D(1) / (4 * D(2).sqrt())
    rows = []
    for i in range(n):
        row = {i: D(i) + D('0.5')}
        for offset in (-3, -1, 1, 3):
            j = i + offset
            if not 0 <= j < n:
                continue
            k = min(i, j) + 1
            if abs(offset) == 1:
                entry = 3 * D(k) * D(k).sqrt()
                if offset < 0:
                    entry = -entry
            else:
                entry = D(k * (k + 1) * (k + 2)).sqrt()
                if offset > 0:
                    entry = -entry
            row[j] = coefficient * entry
        rows.append(row)
    return rows


def scaled_determinant(energy, base):
    """det(T_N-energy)/prod(n+1/2), by seven-band LU without pivoting.

    The scaling is independent of energy and does not change the roots.
    A zero elimination pivot causes an explicit failure: it is not treated
    as an eigenvalue.  No zero pivot occurs in the computation of the article.
    """
    D = Decimal
    rows = [row.copy() for row in base]
    n = len(rows)
    determinant = D(1)
    for k in range(n):
        rows[k][k] -= energy
    for k in range(n):
        pivot = rows[k][k]
        if not pivot:
            raise ArithmeticError('Zero LU pivot: change starting guesses or use pivoting.')
        determinant *= pivot / (D(k) + D('0.5'))
        for i in range(k + 1, min(k + 4, n)):
            if k in rows[i]:
                factor = rows[i].pop(k) / pivot
                for j in range(k + 1, min(k + 4, n)):
                    if j in rows[k]:
                        rows[i][j] = rows[i].get(j, D(0)) - factor * rows[k][j]
    return determinant


def secant_root(base, initial_guess, tolerance=Decimal('1e-55')):
    """Find the root near a coarse seed, with explicit convergence checks."""
    D = Decimal
    a = D(initial_guess) - D('0.0002')
    b = D(initial_guess) + D('0.0002')
    fa, fb = scaled_determinant(a, base), scaled_determinant(b, base)
    for iteration in range(40):
        if not fb:
            return b, iteration
        if fb == fa:
            raise ArithmeticError('Secant denominator vanished.')
        c = b - fb * (b - a) / (fb - fa)
        if abs(c - b) < tolerance:
            if abs(c - D(initial_guess)) > D('0.01'):
                raise ArithmeticError('The root left the neighbourhood of its labelled seed.')
            return c, iteration + 1
        a, b, fa, fb = b, c, fb, scaled_determinant(c, base)
    raise ArithmeticError('Secant iteration did not converge.')


def high_precision_levels(n, digits, initial_guesses):
    """Real levels of the compression of size n to the given number of decimal
    digits, one root near each seed; the roots must increase strictly."""
    with localcontext() as context:
        context.prec = digits
        base = decimal_matrix(n)
        roots = [secant_root(base, guess)[0] for guess in initial_guesses]
        if not all(roots[i] < roots[i + 1] for i in range(len(roots) - 1)):
            raise ArithmeticError('Roots are not distinct and ordered.')
        return roots


def finite_difference_levels(half_width, requested_step, count=3):
    """Fourth-order coordinate discretisation of L at m=1, lambda=3.

    Interior grid in [-half_width,half_width], zero endpoint and exterior
    values.  The latter specify the boundary closure of the five-point
    stencil.  The domain test is supplied separately from the mesh test.
    A fixed starting vector makes the sparse eigensolver deterministic.
    """
    intervals = round(2 * half_width / requested_step)
    step = 2 * half_width / intervals
    x = np.linspace(-half_width, half_width, intervals + 1)[1:-1]
    n = len(x)
    matrix = sp.diags(
        [np.ones(n-2), -16*np.ones(n-1), 30*np.ones(n),
         -16*np.ones(n-1), np.ones(n-2)], [-2, -1, 0, 1, 2], format='csc'
    ) / (24 * step**2)
    matrix = matrix.astype(complex) + sp.diags(0.5*x*x - 0.5j*x**3, format='csc')
    starting_vector = np.exp(-x*x/2).astype(complex)
    values = sla.eigs(matrix, k=count, sigma=0.7, v0=starting_vector,
                      tol=1e-12, maxiter=5000, return_eigenvectors=False)
    return step, values[np.argsort(values.real)]


def check_construction():
    """Check the operator convention and its free limit independently."""
    n = 8
    k = np.arange(1, n+3, dtype=float)
    padded_b = np.diag(np.sqrt(k), -1) - np.diag(np.sqrt(k), 1)
    compressed_cubic = (padded_b @ padded_b @ padded_b)[:n, :n]
    matrix = np.diag(np.arange(n)+0.5) + compressed_cubic/(4*np.sqrt(2))
    np.testing.assert_allclose(ghost_matrix(n), matrix, rtol=3e-15, atol=3e-15)
    np.testing.assert_array_equal(ghost_matrix(n, coupling=0), np.diag(np.arange(n)+0.5))
    small_b = padded_b[:n, :n]
    difference = compressed_cubic - small_b @ small_b @ small_b
    boundary = np.zeros((n, n))
    boundary[-2, -1] = n*np.sqrt(n-1)
    boundary[-1, -2] = -n*np.sqrt(n-1)
    np.testing.assert_allclose(difference, boundary, rtol=3e-15, atol=3e-14)
