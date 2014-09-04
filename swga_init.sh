#!/bin/bash
set -e -x # halt on error, echo commands

cd src/cliquer
make
cd ../..
cp src/cliquer/set_finder .
cp swga/parameters.cfg .






