import os
from setuptools import setup

setup(
    name='PrimerSets',
    version='0.1.0',
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=['swga'],
    package_dir={'swga': 'PrimerSets'},
    url='https://github.com/BrissonEEDS/PrimerSets',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {'console_scripts': ['swga = PrimerSets.scripts.swga:main']},
    test_suite='PrimerSets.test'
)
