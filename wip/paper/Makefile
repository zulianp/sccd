# The article. `make` builds main.pdf from what is committed here and needs no
# Python; `make inputs` regenerates the LaTeX tables from the benchmark CSVs.
LATEXMK ?= latexmk
PYTHON  ?= python3

.PHONY: all inputs check clean distclean supplemental

all: main.pdf supplemental_material.pdf

main.pdf: main.tex preamble.tex arxiv.sty refs.bib $(wildcard sections/*.tex) \
          $(wildcard tikz/*.tex) $(wildcard generated/tables/*.tex)
	$(LATEXMK) -pdf -halt-on-error -interaction=nonstopmode main.tex

# The supplement resolves \cref into the article through xr-hyper, so it needs
# main.aux to exist and be current: main.pdf is a prerequisite, not a courtesy.
supplemental: supplemental_material.pdf
supplemental_material.pdf: supplemental_material.tex preamble.tex main.pdf \
          $(wildcard sections/supp-*.tex) $(wildcard generated/tables/*.tex)
	$(LATEXMK) -pdf -halt-on-error -interaction=nonstopmode supplemental_material.tex

# Rebuild the tables from benchmark/results/*.csv. Only needed after a new sweep;
# the generated fragments are committed so the default target does not run this.
inputs:
	$(PYTHON) scripts/make_inputs.py
	$(PYTHON) scripts/make_counters.py
	$(PYTHON) scripts/make_difficulty.py
	$(PYTHON) scripts/make_symbols.py
	$(PYTHON) scripts/make_competitors.py

# Cross-check the hand-written figures in the prose against the generated tables.
check:
	$(PYTHON) scripts/make_inputs.py --check

clean:
	$(LATEXMK) -c

distclean:
	$(LATEXMK) -C
