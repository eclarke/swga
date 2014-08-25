from distutils.core import setup

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
            'filter-primers = PrimerSets.filter_primers:main',
            'fg-locations = PrimerSets.fg_locations:main',
            'mk-primer-graph = PrimerSets.make_graph:main',
            'find-sets = PrimerSets.find_sets:main',
            'process-sets = PrimerSets.process_sets:main',
        ]
    }
)
