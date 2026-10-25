#!/usr/bin/env python3
# Authors: Davide Batic and Denys Dutykh (Mathematics Department, Khalifa University of
# Science and Technology, Abu Dhabi, UAE).
# Copyright (C) 2026 Davide Batic and Denys Dutykh.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this repository.
"""Finite-matrix checks of the pair weights and of the residue sum.

For the compression T of the ghost generator to dimension d, the spectral
weights of the canonical two-point function are formed from the eigenprojections
of T and the indefinite form of the parity eta, and four identities of the
article are checked: the boundary form of the double commutator [X, [T, X]] of
the compressed coordinate X = B/sqrt(2m); the equality of the sum of the weights
with the expectation of that double commutator in the ground level; the
expression of the combined weight of a complex-conjugate pair through the
energy-weighted form 2 [w, (T - E_0) w]/s_0; and the relation between the norm
of the pair insertion w and the real part of its coefficient. The free case
(lambda=0) gives the sum -1.

Conventions match the article: B has negative upper diagonal, positive lower
one; T is the compression of the full cubic operator, not the cube
of the compressed coordinate.

Run from any directory: python3 code/pair_weights.py
The quantities are printed as JSON.
"""
import json
import numpy as np


def comm(a, b):
    """Commutator of two matrices."""
    return a @ b - b @ a


def model(d, m=1.0, coupling=3.0):
    """Compression T of the generator to dimension d, the coordinate X, the parity
    eta, the compressed B, the matrix R = P B**3 P - (P B P)**3 and the coefficient c."""
    q = np.sqrt(np.arange(1, d + 3, dtype=float))
    b_big = np.diag(q, -1) - np.diag(q, 1)
    b = b_big[:d, :d]
    c = coupling / (6 * (2 * m) ** 1.5)
    r = (b_big @ b_big @ b_big)[:d, :d] - b @ b @ b
    t = m * np.diag(np.arange(d) + 0.5) + c * (b_big @ b_big @ b_big)[:d, :d]
    x = b / np.sqrt(2 * m)
    eta = np.diag((-1.0) ** np.arange(d))
    return t, x, eta, b, r, c


def check(d=5, m=1.0, coupling=3.0):
    """Form the weights, check the four identities and return the quantities."""
    t, x, eta, b, r, c = model(d, m, coupling)
    val, vec = np.linalg.eig(t)
    inv = np.linalg.inv(vec)
    real_indices = np.flatnonzero(abs(val.imag) < 1e-10)
    ground = real_indices[np.argmin(val[real_indices].real)]
    e0 = val[ground].real
    psi = vec[:, ground]
    s0_raw = np.vdot(psi, eta @ psi).real
    psi = psi / np.sqrt(abs(s0_raw))
    s0 = np.sign(s0_raw)
    bracket = lambda a, b: np.vdot(a, eta @ b)
    insertion = x @ psi
    projections = [np.outer(vec[:, n], inv[n, :]) for n in range(d)]
    weights = np.array([2 * (val[n] - e0) * bracket(insertion, p @ insertion) / s0 for n, p in enumerate(projections)])
    total = np.sum(weights)
    double = comm(x, comm(t, x))
    direct = bracket(psi, double @ psi) / s0
    boundary = -np.eye(d)
    boundary[-1, -1] += d
    boundary += c / (2 * m) * comm(b, comm(r, b))
    pair = int(np.argmax(val.imag))
    mate = int(np.argmin(abs(val - np.conj(val[pair]))))
    w = (projections[pair] + projections[mate]) @ insertion
    cp = bracket(w, projections[pair] @ w) / s0
    norm = bracket(w, w).real
    energy_form = 2 * bracket(w, (t - e0 * np.eye(d)) @ w) / s0
    assert np.max(abs(double - boundary)) < 2e-12
    assert abs(total - direct) < 2e-12
    assert abs(energy_form - weights[pair] - weights[mate]) < 2e-12
    assert abs(norm / s0 - 2 * cp.real) < 2e-12
    return {"dimension": d, "m": m, "lambda": coupling,
            "E0": e0, "pair_energy": [val[pair].real, val[pair].imag],
            "pair_insertion_norm": norm,
            "combined_pair_weight": (weights[pair] + weights[mate]).real,
            "spectral_weight_sum": total.real,
            "double_commutator_expectation": direct.real,
            "boundary_identity_max_error": float(np.max(abs(double - boundary)))}


if __name__ == "__main__":
    result = check()
    t, x, eta, *_ = model(5, 1.0, 0.0)
    psi0 = np.eye(5)[:, 0]
    free_sum = np.vdot(psi0, eta @ comm(x, comm(t, x)) @ psi0).real
    assert abs(free_sum + 1) < 1e-14
    result["free_sum"] = free_sum
    print(json.dumps(result, indent=2))
