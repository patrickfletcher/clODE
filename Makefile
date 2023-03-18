
PYFILES=$(shell find clode/python -name "*.py")
PYTESTFILES=$(shell find test -name "*.py")

venv:
	python -m venv venv

install:
	. venv/bin/activate && pip install --upgrade pip && \
	python -m pip install -r requirements.txt

install_clode:
	. venv/bin/activate && python -m pip install .

all: install install_clode
	. venv/bin/activate && ls

format:
	. venv/bin/activate && \
		isort $(PYFILES) $(PYTESTFILES) && \
		black $(PYFILES) $(PYTESTFILES)

test: install install_clode
	echo $(PYTESTFILES)
	. venv/bin/activate && python -m pytest $(PYTESTFILES)

run: install
	. venv/bin/activate && PYTHONPATH=$(PYTHONPATH) python main.py

lint: install
	. venv/bin/activate && vulture $(PYFILES) $(PYTESTFILES) && \
	python -m pylint $(PYFILES) $(PYTESTFILES) && \
	mypy $(PYFILES) $(PYTESTFILES)