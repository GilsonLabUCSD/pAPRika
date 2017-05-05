#!/bin/sh
set -ex

echo "Simplest possible testing regimen"
cd test
python test_alignment.py
python test_solvation.py
