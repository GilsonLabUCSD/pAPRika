#!/bin/sh
set -ex

echo "Run all tests except slow ones"
cd test && pytest -v -x -m 'not slow' --cov=paprika .
