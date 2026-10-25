# Makefile of the numerical material accompanying the article
# "Spectrum and positive metrics of the cubic ghost oscillator"
# Authors: Davide Batic and Denys Dutykh (Mathematics Department, Khalifa University of
# Science and Technology, Abu Dhabi, UAE)
# Copyright (C) 2026 Davide Batic and Denys Dutykh.
# Distributed under the GNU Lesser General Public License, version 2.1; see LICENSE.
#
# make          the tables, the figure, then the pair-weight check
# make spectra  the tables of data/ and the figure of figures/
# make pairs    the finite-matrix check of the pair weights and of the sum rule
# make clean    remove the Python cache; the tables and the figure stay

PYTHON ?= python3

.PHONY: all spectra pairs clean

all: spectra pairs

spectra:
	$(PYTHON) code/truncation_spectra.py

pairs:
	$(PYTHON) code/pair_weights.py > data/pair_weights.json

clean:
	rm -rf code/__pycache__
