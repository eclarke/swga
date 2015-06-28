from setuptools import setup, find_packages

with open("VERSION") as version_file:
    version = version_file.read().strip()

setup(
    name='swga',
    version=version,
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=find_packages(),
    install_requires=['pyfaidx',
                      'click',
                      'pyyaml',
                      'peewee',
                      'data_hacks',
                      'pytest'],
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
    entry_points = {'console_scripts': ['swga = swga.commands.swga_entry:main']}
)
