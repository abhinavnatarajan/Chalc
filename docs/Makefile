# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

stubs:
	@echo "Generating stubs for chalc.chromatic"
	@python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir ../src
	@echo "Generating stubs for chalc.filtration"
	@python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir ../src

.PHONY: help Makefile stubs

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile stubs
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
