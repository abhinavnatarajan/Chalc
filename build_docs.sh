#!/bin/bash

python docs/stubgen.py chalc --output-dir docs/stubs
shopt -s extglob
cp src/chalc/!(__init__)*.py docs/stubs/chalc/ -rf

sphinx-build -E -b html docs docs/_build