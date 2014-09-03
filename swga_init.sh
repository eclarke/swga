#!/bin/bash
set -e -x # halt on error, echo commands

cd lib
make
cd ..
cp lib/set_finder .
cp PrimerSets/parameters.cfg ..






