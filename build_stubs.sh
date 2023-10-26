#!/bin/bash

python -m pybind11_stubgen chalc.chromatic --numpy-array-wrap-with-annotated --output-dir src
python -m pybind11_stubgen chalc.filtration --numpy-array-wrap-with-annotated --output-dir src
