#!/bin/bash

python docs/stubgen.py chalc --output-dir docs/stubs
cp src/chalc/sixpack/* docs/stubs/chalc/sixpack -rf

sphinx-build -E -b html docs docs/_build