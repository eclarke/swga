import os
from setuptools import setup, find_packages

setup(
    name='swga',
    version='0.2.1',
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=find_packages(exclude=['cliquer', 'docs', 'tests']),
    install_requires=['pyfasta'],
    url='https://github.com/BrissonEEDS/PrimerSets',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {'console_scripts': ['swga = swga.commands.swga_entry:main']},
)
