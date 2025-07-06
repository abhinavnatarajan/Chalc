# This makefile is used to generate type stubs for chalc.chromatic and chalc.filtration.
# Put it first so that "make" without argument is like "make stubs".
stubs:
	@echo "Generating stubs for chalc.chromatic"
	@python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ./src
	@echo "Generating stubs for chalc.filtration"
	@python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ./src

.PHONY: Makefile stubs all clean test

clean:
	rm src/chalc/chromatic.pyi src/chalc/filtration.pyi
