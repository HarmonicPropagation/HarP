# Reusable Python environment Makefile

VENV := .venv
PYTHON := $(VENV)/bin/python3
PIP := $(VENV)/bin/pip
PIP_COMPILE := $(VENV)/bin/pip-compile
PIP_SYNC := $(VENV)/bin/pip-sync

# 1. Global exports & OS-aware variables
# Ensure local modules are discoverable for all test/typecheck targets
export PYTHONPATH := ./src

# macOS (Darwin) uses 'open', Linux uses 'xdg-open'
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	OPEN_CMD := open
else
	OPEN_CMD := xdg-open
endif

# 2. Phony targets
.PHONY: help init compile update clean compile-all dev docs view-docs pdf-docs typecheck test benchmark clean-tests

help:
	@echo "Available commands:"
	@echo "  make init        - Create .venv and install base dependencies"
	@echo "  make dev         - Create/sync .venv with development and testing dependencies"
	@echo "  make compile     - Generate requirements.txt from requirements.in"
	@echo "  make compile-all - Generate both base and dev requirements"
	@echo "  make update      - Sync .venv to match requirements.txt"
	@echo "  make test        - Run standard correctness and accuracy tests"
	@echo "  make benchmark   - Run performance and complexity benchmarks"
	@echo "  make typecheck   - Run static type checking via mypy"
	@echo "  make docs        - Build HTML documentation via Sphinx"
	@echo "  make view-docs   - Build and open the HTML documentation"
	@echo "  make clean       - Remove .venv, caches, and test artifacts"

# ---------------------------------------------------------
# Base environment & compilation
# ---------------------------------------------------------

$(VENV)/bin/activate:
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip pip-tools

requirements.txt: requirements.in | $(VENV)/bin/activate
	$(PIP_COMPILE) --strip-extras requirements.in

init: requirements.txt | $(VENV)/bin/activate
	$(PIP_SYNC) requirements.txt
	@echo "\nEnvironment ready. Run 'source $(VENV)/bin/activate' to enter."

compile: requirements.txt
update: init

# ---------------------------------------------------------
# Development & analysis environment
# ---------------------------------------------------------

requirements-dev.txt: requirements-dev.in requirements.txt | $(VENV)/bin/activate
	$(PIP_COMPILE) --strip-extras requirements-dev.in

compile-all: requirements.txt requirements-dev.txt

dev: requirements-dev.txt | $(VENV)/bin/activate
	$(PIP_SYNC) requirements.txt requirements-dev.txt
	@echo "\nDev environment ready. Notebook dependencies installed."

# ---------------------------------------------------------
# Documentation
# ---------------------------------------------------------

docs: dev
	$(MAKE) -C doc html SPHINXBUILD=$(abspath $(VENV)/bin/sphinx-build)

view-docs: docs
	$(OPEN_CMD) doc/build/html/index.html

pdf-docs: dev
	$(MAKE) -C doc latexpdf SPHINXBUILD=$(abspath $(VENV)/bin/sphinx-build)

# ---------------------------------------------------------
# Testing, benchmarking, and type checking
# ---------------------------------------------------------

typecheck: dev
	$(VENV)/bin/mypy src/ tests/ --ignore-missing-imports --explicit-package-bases --strict

test: dev
	$(VENV)/bin/pytest -v tests/

benchmark: dev
	$(VENV)/bin/pytest --benchmark-only tests/test_performance.py

# ---------------------------------------------------------
# Cleanup
# ---------------------------------------------------------

clean: clean-tests
	rm -rf $(VENV)
	find . -type d -name "__pycache__" -exec rm -rf {} +

clean-tests:
	rm -rf .pytest_cache/
	rm -rf .benchmarks/
	find tests/ -type f -name "*.h5" -delete
