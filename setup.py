from setuptools import setup, find_packages, Command
from distutils.command.build import build
from subprocess import call
import os, sys

import swga

basepath = os.path.dirname(os.path.abspath(__file__))
cl_path = os.path.join(basepath, "swga/contrib/cliquer")

class build_swga(build):
    def run(self):
        def compile():
            call(['make'], cwd=cl_path)
        build.run(self)
        self.execute(compile, [], "Compiling needed binaries...")

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
    packages=['swga'],
    install_requires=['pyfaidx',
                      'click'],
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],
    package_data={'swga': ['data/*']},
    url='https://github.com/eclarke/swga',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {'console_scripts': ['swga = swga.commands.swga_entry:main']},
    cmdclass={
        'build': build_swga,
        'test': test_swga
    }
)
