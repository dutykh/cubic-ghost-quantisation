#!/usr/bin/env python3
# Authors: Davide Batic and Denys Dutykh (Mathematics Department, Khalifa University of
# Science and Technology, Abu Dhabi, UAE).
# Copyright (C) 2026 Davide Batic and Denys Dutykh.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this repository.
"""Spectra of the finite compressions of the cubic ghost generator.

Computes the tables behind the figure of finite-truncation spectra and draws
the figure: the eigenvalues of the compressions of sizes 40 and 80 with their
condition numbers and residuals; the 31 lowest real levels from the
compressions of sizes 300 and 360 in 100- and 110-digit arithmetic, and a check
of the highest of them at 140 digits; and a fourth-order finite-difference
calculation in coordinate space. The starting values of the levels are read from
data/starting_guesses.csv; the tables are written to data/ and the figure to
figures/. No random numbers are used.

Run from any directory: python3 code/truncation_spectra.py
The option --skip-fd omits the coordinate-space calculation.
"""
import argparse
import csv
from decimal import Decimal, localcontext
import logging
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from ghost_compressions import (check_construction, spectrum_and_condition_numbers,
                                high_precision_levels, finite_difference_levels)

# The font subsetter warns about the creation dates of the bundled fonts; the warning is harmless.
logging.getLogger('fontTools').setLevel(logging.ERROR)

ROOT = Path(__file__).resolve().parents[1]


def write_csv(path, header, rows):
    """Write a table with its header row."""
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def plot_figure(values, reference, output):
    """Draw the spectra of the two compressions against the reference levels."""
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm',
                         'font.size': 8, 'axes.labelsize': 9,
                         'axes.titlesize': 9, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'svg.hashsalt': 'ghost-truncation-spectrum'})
    fig, axes = plt.subplots(2, 1, figsize=(6.05, 4.3), sharex=True)
    colors = ['#8E2A3F', '#17365D']
    markers = ['s', 'o']
    for panel, (ax, dimension) in enumerate(zip(axes, [40, 80])):
        ax.set_yscale('symlog', linthresh=1, linscale=0.6, base=10)
        ax.axhline(0, color='#bfc3c6', linewidth=0.45, zorder=1)
        ax.plot(reference, np.zeros_like(reference), '|', color='#5C6670',
                markersize=12, markeredgewidth=0.65, zorder=2)
        w = values[dimension]
        ax.plot(w.real, w.imag, linestyle='none', marker=markers[panel],
                markerfacecolor='white', markeredgecolor=colors[panel],
                markeredgewidth=0.8, markersize=4, zorder=3)
        ax.set_xlim(0, 77)
        ax.set_ylim(-1600, 1600)
        ax.set_yticks([-1000, -100, -10, -1, 0, 1, 10, 100, 1000])
        ax.set_yticklabels([r'$-10^3$', r'$-10^2$', r'$-10^1$', r'$-1$', '$0$',
                            '$1$', r'$10^1$', r'$10^2$', r'$10^3$'])
        ax.set_ylabel(r'$\operatorname{Im}E$')
        ax.set_title(f'({chr(97+panel)}) $N={dimension}$', loc='left', pad=3)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(width=0.6, length=3)
        for spine in ax.spines.values():
            spine.set_linewidth(0.6)
    axes[-1].set_xlabel(r'$\operatorname{Re}E$')
    handles = [Line2D([], [], color='#5C6670', marker='|', markersize=7,
                      linestyle='none', label=r'reference levels $E_n$'),
               Line2D([], [], color=colors[0], marker='s', markersize=4,
                      markerfacecolor='white', linestyle='none', label=r'$N=40$'),
               Line2D([], [], color=colors[1], marker='o', markersize=4,
                      markerfacecolor='white', linestyle='none', label=r'$N=80$')]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(0.55, 0.995), columnspacing=1.6, handletextpad=0.45)
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.09, top=0.89, hspace=0.22)
    metadata = {'Title': 'Finite compressions of the cubic ghost generator',
                'Author': 'Davide Batic and Denys Dutykh',
                'Subject': 'm=1, lambda=3; P_J B^3 P_J; N=J+1',
                'Creator': None, 'Producer': None,
                'CreationDate': None, 'ModDate': None}
    fig.savefig(output / 'ghost-truncation-spectrum.pdf', metadata=metadata)
    fig.savefig(output / 'ghost-truncation-spectrum.svg',
                metadata={'Title': metadata['Title'], 'Date': None, 'Creator': None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-fd', action='store_true', help='Omit the coordinate-space check.')
    args = parser.parse_args()
    data, figures = ROOT/'data', ROOT/'figures'
    data.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    check_construction()
    values = {}
    for n in [40, 80]:
        w, condition, residual = spectrum_and_condition_numbers(n)
        values[n] = w
        write_csv(data/f'spectrum_N{n}.csv',
                  ['index', 'real_energy', 'imag_energy', 'eigenvalue_condition_number',
                   'relative_right_residual'],
                  [[j, format(z.real, '.17g'), format(z.imag, '.17g'),
                    format(condition[j], '.17g'), format(residual[j], '.17g')]
                   for j, z in enumerate(w)])
        print(f'N={n}: maximum relative residual {max(residual):.3e}', flush=True)

    with (data/'starting_guesses.csv').open() as handle:
        guesses = [row['initial_energy_guess'] for row in csv.DictReader(handle)]
    roots = {}
    for n, digits in [(300, 100), (360, 110)]:
        roots[n] = high_precision_levels(n, digits, guesses)
        print(f'N={n}: {len(roots[n])} real roots at {digits} digits.', flush=True)
    with localcontext() as context:
        context.prec = 120
        differences = [abs(a-b) for a, b in zip(roots[300], roots[360])]
        write_csv(data/'reference_levels.csv',
                  ['n', 'energy_N300_100digit', 'energy_N360_110digit', 'absolute_difference'],
                  [[n, str(a), str(b), str(differences[n])]
                   for n, (a,b) in enumerate(zip(roots[300], roots[360]))])
        max_difference = max(differences)
        assert max_difference < Decimal('1e-8')
        print('Maximum N=300 versus N=360 difference:',
              format(max_difference, '.6e'), flush=True)

    # Additional arithmetic check at the highest displayed level.
    arithmetic = []
    with localcontext() as context:
        context.prec = 155
        for n in [300, 360]:
            check = high_precision_levels(n, 140, [guesses[30]])[0]
            delta = abs(check-roots[n][30])
            arithmetic.append([n, 30, str(check), str(delta)])
    write_csv(data/'arithmetic_check.csv',
              ['dimension', 'level', 'energy_140digit', 'difference_from_figure_precision'], arithmetic)

    plot_figure(values, np.array([float(x) for x in roots[360]]), figures)
    if not args.skip_fd:
        fd_rows = []
        # Mesh refinement at L=10, then a domain-size comparison at h=.02.
        configurations = [(10.,h,'mesh') for h in [.04,.02,.01,.005]]
        configurations += [(L,.02,'domain') for L in [8.,12.]]
        for half_width, requested_step, test in configurations:
            step, w = finite_difference_levels(half_width, requested_step)
            for level, z in enumerate(w):
                error = abs(z-complex(float(roots[360][level]), 0))
                fd_rows.append([test, half_width, step, level,
                                format(z.real,'.17g'), format(z.imag,'.17g'),
                                format(error,'.17g')])
            print(f'FD half-width={half_width:g}, h={step:g}: E0={w[0]}', flush=True)
        write_csv(data/'coordinate_check.csv',
                  ['test', 'box_half_width', 'step', 'level', 'real_energy', 'imag_energy',
                   'absolute_error_vs_N360'], fd_rows)


if __name__ == '__main__':
    main()
