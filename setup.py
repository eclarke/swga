import os
import re
import platform
import subprocess
from setuptools import setup, find_packages, Command, dist
from distutils.command.build import build as _build
from distutils.command.build_py import build_py as _build_py
from distutils.spawn import find_executable
from setuptools.command.bdist_egg import bdist_egg as _bdist_egg

version = "0.4.0"

class BinaryDistribution(dist.Distribution):
    def is_pure(self):
        return False



class bdist_egg(_bdist_egg):
    def run(self):
        self.run_command('build_py')
        _bdist_egg.run(self)


class build_py(_build_py):

    def run(self):
        """
        Builds `DSK` and `cliquer` in the contrib/ directory and moves them 
        to the right place before installing `swga`.
        """
        def find_prog(target='g++'):
            '''Finds all items on the path with the target in their name.'''
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip("")
                if not os.path.isdir(path):
                    continue
                for item in os.listdir(path):
                    if target in item:
                        yield item

        def get_version(target, string):
            '''Finds a version number (i.e. gcc-4.9 ==> 4.9)'''
            match = re.match(r'{}-([0-9.]+)'.format(re.escape(target)), string)
            if match:
                return float(match.group(1))

        # Defaults (e.g. for Linux)
        CC = find_executable('g++')
        omp = 1
        osx = 0
        # Check if we're on a Mac and adjust compiler accordingly
        if 'Darwin' in platform.system() and not os.environ.get("CC"):
            osx = 1
            gcc_version = subprocess.check_output(
                [CC, '--version'], 
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
                        "of clang. Override by setting CC=/path/to/gcc before "
                        "this command.\n".format(real_CC)
                    )
                    CC = real_CC
                else:
                    print(
                        "\nFound clang but no GNU g++ compiler: Disabling "
                        "multiprocessing for find_sets. If GNU g++ is "
                        "installed, set it with CC=/path/to/gcc before this "
                        "command."
                    )
                    omp=0

        cliquer_make_cmd = (['make'])
        dsk_make_cmd = (
            'make',
            'CC='+CC,
            'omp='+str(omp),
            'osx='+str(osx)
        )

        print(" ".join(cliquer_make_cmd))
        if not self.dry_run:
            subprocess.check_call(cliquer_make_cmd, cwd="contrib/cliquer")
        print(" ".join(dsk_make_cmd))
        if not self.dry_run:
            subprocess.check_call(dsk_make_cmd, cwd="contrib/dsk")

        if not self.dry_run:            
            try:
                subprocess.check_call(["mkdir", "-p", "swga/bin"])
                print("Moving set_finder")
                os.rename("contrib/cliquer/set_finder", "swga/bin/set_finder")
                print("Moving dsk")
                os.rename("contrib/dsk/dsk", "swga/bin/dsk")
            except OSError:
                raise

        _build_py.run(self)

setup(
    name='swga',
    version=version,
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=find_packages(),
    install_requires=[
        'pyfaidx',
        'click',
        'pyyaml',
        'peewee',
        'melt',
        'pytest',
        'argutils',
        'semantic_version'
    ],
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],
    package_data={
        'swga': ['commands/specfiles/*', 'bin/set_finder', 'bin/dsk']
    },
    data_files=[ 
        ('bin', ['swga/bin/set_finder', 'swga/bin/dsk'])
    ],
    eager_resources=["bin/set_finder", "bin/dsk"],
    include_package_data=True,
    distclass=BinaryDistribution,
    url='https://github.com/eclarke/swga',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {'console_scripts': ['swga = swga.commands.swga_entry:main']},
    cmdclass={
        'bdist_egg': bdist_egg,
        'build_py': build_py
    }
)
