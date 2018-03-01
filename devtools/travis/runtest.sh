#!/bin/sh
set -ex

echo "Simplest possible testing regimen"
cd test && pytest -v -x --cov=paprika .
