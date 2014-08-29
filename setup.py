from setuptools import setup

setup(
    name='PrimerSets',
    version='0.1.0',
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=['PrimerSets'],
    scripts=['lib/set_finder', 'utils/fasta_flattener.sh'],
    url='https://github.com/BrissonEEDS/PrimerSets',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {
        'console_scripts': [
            'filter_primers = PrimerSets.scripts.filter_primers:main',
            'find_fg_locations = PrimerSets.scripts.fg_locations:main',
            'mk_primer_graph = PrimerSets.scripts.make_graph:main',
            'find_sets = PrimerSets.scripts.find_sets:main',
            'process_sets = PrimerSets.scripts.process_sets:main',
        ]
    },
    test_suite='PrimerSets.test'
)
