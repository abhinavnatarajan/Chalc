#!/bin/bash

python -m pybind11_stubgen chalc.chromatic --numpy-array-use-type-var --output-dir src
python -m pybind11_stubgen chalc.filtration --numpy-array-use-type-var --output-dir src
sphinx-build -E -b html docs docs/_build
