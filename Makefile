# This makefile is used to generate type stubs for chalc.chromatic and chalc.filtration.
# Put it first so that "make" without argument is like "make stubs".
install: dev stubs lock

dev: # install editable project + dev/docs/test dependencies
	uv sync --verbose --group dev --group docs --group test --no-progress

upgrade: # update project dependencies
	uv sync -U --group dev --group docs --group test --verbose --no-progress

lock: # update lockfiles
	uv lock
	uv export --locked --format pylock.toml --group dev --group docs --group test -o pylock.toml --quiet

wheel: stubs lock # build the project wheel
	uv build --no-progress --verbose --wheel

sdist: # build source distribution
	uv build --no-progress --verbose --sdist

tests: # run tests without requiring GUDHI
	uv run --group test pytest

tests-gudhi: # include the optional GUDHI comparisons
	uv run --group test-gudhi pytest

stubs: dev # stubs
	@echo 'Generating stubs for chalc.chromatic'
	uv run python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ./src
	@echo 'Generating stubs for chalc.filtration'
	uv run python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ./src

docs: stubs
	$(MAKE) -C docs html

all: install wheel sdist docs

clean:
	rm src/chalc/chromatic.pyi src/chalc/filtration.pyi
	$(MAKE) -C docs clean

.PHONY: install dev upgrade lock wheel sdist tests tests-gudhi stubs clean
