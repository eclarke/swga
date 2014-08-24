from distutils.core import setup

setup(
    name='PrimerSets',
    version='0.1.0',
    author='Erik Clarke',
    author_email='ecl@mail.med.upenn.edu',
    packages=['swga'],
    scripts=['lib/set_finder', 'utils/fasta_flattener.sh'],
    url='https://github.com/BrissonEEDS/PrimerSets',
    license='LICENSE.txt',
    description='Pipeline to select compatible primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points = {
        'console_scripts': [
            'filter_primers = swga.filter_primers:main',
            'fg_locations = swga.fg_locations:main',
            'mk_primer_graph = swga.make_graph:main',
            'find_sets = swga.find_sets:main',
            'process_sets = swga.process_sets:main',
        ]
    }
)
    
    
