# -*- coding: utf-8 -*-
import os
from swga.commands2 import Command
import swga.primers as primers

def main(argv, cfg_file):
    cmd = Command('count', cfg_file=cfg_file)
    cmd.parse_args(argv)
    count_kmers(**cmd.args)

    
def count_kmers(fg_genome_fp,
                bg_genome_fp, 
                min_size, 
                max_size, 
                threshold, 
                output_dir):
    assert os.path.isfile(fg_genome_fp)
    assert os.path.isfile(bg_genome_fp)

    for i in xrange(min_size, max_size + 1):
        primers.count_kmers(i, fg_genome_fp, output_dir, threshold)
        primers.count_kmers(i, bg_genome_fp, output_dir, threshold)

