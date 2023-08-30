#!/bin/bash

python stubgen.py chalc

sphinx-build -E -b html docs/source docs/build