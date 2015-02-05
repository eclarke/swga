from setuptools import setup, Command, find_packages
from subprocess import call
import sys

import swga

class test_swga(Command):
    user_options = []
    def initialize_options(self): pass
    def finalize_options(self): pass

    def run(self):
        errno = call([sys.executable, 'test_swga.py'])
        raise SystemExit(errno)

setup(
    name='swga',
    version=swga.__version__,
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=find_packages(),
    install_requires=['pyfaidx',
                      'click',
                      'pyyaml'],
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],
    package_data={'swga': ['data/*',
                           'bin/*']},
    url='https://github.com/eclarke/swga',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {'console_scripts': ['swga = swga.commands.swga_entry:main']},
    cmdclass={
        'test': test_swga
    }
)
