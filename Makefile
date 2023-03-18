
PYFILES=$(shell find clode/python -name "*.py")
PYTESTFILES=$(shell find test -name "*.py")

venv:
	python -m venv venv

install:
	pip install --upgrade pip && \
		python -m pip install -r requirements.txt

install_clode:
	python -m pip install .

format:
	isort $(PYFILES) $(PYTESTFILES) && \
			black $(PYFILES) $(PYTESTFILES)

test: install install_clode
	python -m pytest $(PYTESTFILES)

run: install
	. venv/bin/activate && PYTHONPATH=$(PYTHONPATH) python main.py

lint: install
	vulture $(PYFILES) $(PYTESTFILES) && \
		python -m pylint $(PYFILES) $(PYTESTFILES) && \
		mypy $(PYFILES) $(PYTESTFILES)