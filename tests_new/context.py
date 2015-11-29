"""
Allows the tests to import the source code within this directory without
running `setup.py develop` or `pip install --editable`.
From Kenneth Reitz (http://www.kennethreitz.org/repository-structure-and-python/)
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import swga