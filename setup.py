import os
import sys
import shutil
import stat
import re
import platform
import subprocess
from setuptools import setup, find_packages, dist
from distutils.command.build_py import build_py as _build_py
from distutils.spawn import find_executable
from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
from setuptools.command.develop import develop as _develop

version = "0.4.4"

# Test to ensure Python < 3 (if pip doesn't check for it later)
if not sys.version_info[0] == 2:
    sys.exit("Python 3 is not supported (sorry)")

class develop(_develop):

    def run(self):
        self.run_command('build_py')
        _develop.run(self)

class bdist_egg(_bdist_egg):

    def run(self):
        self.run_command('build_py')
        _bdist_egg.run(self)


class build_py(_build_py):

    def run(self):
        """
        Builds `DSK` and `cliquer` in the ext/ directory and moves them 
        to the right place before installing `swga`.
        """
        def find_prog(target='g++'):
            """Finds all items on the path with the target in their name."""
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip("")
                if not os.path.isdir(path):
                    continue
                for item in os.listdir(path):
                    if target in item:
                        yield item

        def get_version(target, string):
            """Finds a version number (i.e. gcc-4.9 ==> 4.9)"""
            match = re.match(r'{}-([0-9.]+)'.format(re.escape(target)), string)
            if match:
                return float(match.group(1))

        def make_executable(fp):
            st = os.stat(fp)
            os.chmod(fp, st.st_mode | stat.S_IEXEC)

        # Defaults (e.g. for Linux)
        gcc = find_executable('g++')
        omp = 1
        osx = 0
        # Check if we're on a Mac and adjust compiler accordingly
        if 'Darwin' in platform.system() and not os.environ.get("CC"):
            osx = 1
            gcc_version = subprocess.check_output(
                [gcc, '--version'],
                stderr=subprocess.STDOUT
            )
            if 'LLVM' in gcc_version:
                # We'd rather not use clang, so search the path for the most
                # recent version of g++ we can find
                bestv = (0, '')
                for hit in find_prog("g++"):
                    ver = get_version("g++", hit)
                    if ver >= bestv[0]:
                        bestv = (ver, hit)
                real_CC = find_executable(bestv[1])
                if real_CC:
                    print(
                        "\nNote: Found real g++ living in `{}`, using it instead "
                        "of clang.\nOverride by setting CC=/path/to/gcc before "
                        "this command.\n".format(real_CC)
                    )
                    gcc = real_CC
                else:
                    print(
                        "\nFound clang but no GNU g++ compiler: Disabling "
                        "multiprocessing for find_sets. If GNU g++ is "
                        "installed, set it with CC=/path/to/gcc before this "
                        "command."
                    )
                    omp = 0

        cliquer_make_cmd = (['make'])
        dsk_make_cmd = (
            'make',
            'CC=' + gcc,
            'omp=' + str(omp),
            'osx=' + str(osx)
        )

        if not self.dry_run:
            subprocess.check_call(cliquer_make_cmd, cwd="ext/cliquer")
            subprocess.check_call(dsk_make_cmd, cwd="ext/dsk")
            make_executable("ext/cliquer/set_finder")
            make_executable("ext/dsk/dsk")
            if not os.path.exists('swga/bin'):
                os.mkdir('swga/bin')
            shutil.copyfile("ext/cliquer/set_finder", "swga/bin/set_finder")
            shutil.copyfile("ext/dsk/dsk", "swga/bin/dsk")
            make_executable("swga/bin/set_finder")
            make_executable("swga/bin/dsk")

        _build_py.run(self)

setup(
    name='swga',
    version=version,
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=find_packages(),
    install_requires=[
        'pyfaidx>0.4.5.2',
        'click',
        'pyyaml',
        'peewee>=2.7.3,<3.0',
        'melt',
        'pytest',
        'argutils',
        'semantic_version'
    ],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],
    package_data={
        'swga': ['commands/specfiles/*', 'bin/*']
    },
    data_files=[
        ('bin', ['ext/dsk/dsk', 'ext/cliquer/set_finder'])
    ],
    include_package_data=True,
    url='https://github.com/eclarke/swga',
    license='LICENSE.txt',
    description='Pipeline to select primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points={'console_scripts': ['swga = swga.main:main']},
    cmdclass={
        'bdist_egg': bdist_egg,
        'build_py': build_py,
        'develop': develop
    }
)
