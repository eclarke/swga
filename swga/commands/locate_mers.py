import os
import sys
import multiprocessing
import swga
import swga.primers
import swga.genome
from swga.commands import Command

def main(argv, cfg_file):
    cmd = Command('locate-mers', cfg_file=cfg_file)
    cmd.parse_args(argv)
    locate_primers(**cmd.args)

def locate_primers(input,
                   output,
                   genome,
                   ncores,
                   passthrough):

    # Check to make sure genome file is valid
    if not os.path.isfile(genome):
        swga.swga_error('Error: Genome specified by %s does not exist.' %
                         genome)
    if not swga.genome.check_if_flattened(genome):
        swga.swga_error('Error: Genome does not appear to be flattened: use '
                         '`swga flatten` first.')

    primers = swga.primers.read_primer_file(input, passthrough)


    # Find locations using multiple processes
    locations = swga.genome.mp_find_primer_locations(primers, genome, ncores)

    # Save to gzipped pickled file (optimized for large numbers of sites)
    sys.stderr.write("Saving primer locations...\n")
    swga.genome.save_locations(locations, output)
