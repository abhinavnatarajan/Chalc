# This makefile is used to generate type stubs for chalc.chromatic and chalc.filtration.
# Put it first so that "make" without argument is like "make stubs".
install: dev stubs lock

dev: # install editable project + dev/test dependencies
	uv sync --verbose --all-groups --no-progress

upgrade: # update project dependencies
	uv sync -U --all-groups --verbose --no-progress

lock: # update lockfiles
	uv lock
	uv export --format pylock.toml --all-groups -o pylock.toml --quiet

wheel: stubs lock # build the project wheel
	uv build --no-progress --verbose --wheel

sdist: # build source distribution
	uv build --no-progress --verbose --sdist

tests: dev # run the test suite
	uv run pytest

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

.PHONY: install dev upgrade lock wheel sdist tests stubs clean
