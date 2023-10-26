#!/bin/bash

# python docs/stubgen.py chalc --output-dir docs/stubs
source build_stubs.sh
shopt -s extglob
mkdir -p docs/stubs/chalc/
cp -rf src/chalc/!(__init__)*.py?(i) docs/stubs/chalc/
for f in docs/stubs/chalc/*.pyi; do mv -- "$f" "${f%.pyi}.py"; done
sphinx-build -E -b html docs docs/_build
