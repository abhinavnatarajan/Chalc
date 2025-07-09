# This makefile is used to generate type stubs for chalc.chromatic and chalc.filtration.
# Put it first so that "make" without argument is like "make stubs".
stubs:
	@echo "Generating stubs for chalc.chromatic"
	@python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ./src
	@echo "Generating stubs for chalc.filtration"
	@python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ./src

.PHONY: stubs all clean test docs

docs:
	cd docs && make html

all: stubs docs

clean:
	rm src/chalc/chromatic.pyi src/chalc/filtration.pyi
	cd docs && make clean
