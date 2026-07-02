PYTHON ?= python3
CHROMO = PYTHONPATH=src $(PYTHON) -m chromosort.cli

.PHONY: help test smoke diff-check lint agent-check docs-check package-check

help:
	@printf '%s\n' \
		'ChromoSort development commands:' \
		'  make test          Run the synthetic unit test suite.' \
		'  make smoke         Run CLI help smoke checks.' \
		'  make diff-check    Check for whitespace errors in the git diff.' \
		'  make lint          Run conservative Ruff checks.' \
		'  make agent-check   Run diff-check, lint, smoke, test, and package-check.' \
		'  make docs-check    Build the Jekyll docs site if bundle is available.' \
		'  make package-check Build sdist/wheel and verify import metadata.'

test:
	$(PYTHON) -m unittest discover -s tests -v

smoke:
	$(CHROMO) --help
	$(CHROMO) sort --help
	$(CHROMO) clean --help
	$(CHROMO) eval all --help
	$(CHROMO) fix --help
	$(CHROMO) cut --help
	$(CHROMO) manual --help
	$(CHROMO) gafprep --help
	$(CHROMO) plot --help
	$(CHROMO) graph-map --help
	$(CHROMO) scaffold --help
	$(CHROMO) gapfill --help

diff-check:
	git diff --check

lint:
	ruff check .

agent-check: diff-check lint smoke test package-check

docs-check:
	bundle exec jekyll build --source docs --destination docs/_site

package-check:
	$(PYTHON) -m build --sdist --wheel --no-isolation
	PYTHONPATH=src $(PYTHON) -c "import chromosort; print(chromosort.__version__)"
