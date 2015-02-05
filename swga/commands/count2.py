# -*- coding: utf-8 -*-
import os
import swga
import swga.primers as primers
from swga.primers import Primer
from swga.commands2 import Command
from swga.melting import Tm
from swga.clint.textui import progress
import click

def main(argv, cfg_file):
    cmd = Command('count', cfg_file=cfg_file)
    cmd.parse_args(argv)
    count_kmers(**cmd.args)

    
def count_kmers(fg_genome_fp,
                bg_genome_fp, 
                min_size, 
                max_size, 
                threshold, 
                primer_db,
                exclude_fp,
                exclude_threshold):
    assert os.path.isfile(fg_genome_fp)
    assert os.path.isfile(bg_genome_fp)

    if os.path.isfile(primer_db):
        swga.warn("Existing database found at %s" % os.path.abspath(primer_db))
        swga.warn("Re-counting primers will reset the entire database!")
        click.confirm("Are you sure you want to proceed?", abort=True)

    primers.init_db(primer_db, create_if_missing=True)
    primers.create_tables()
    
    output_dir = ".swga_tmp"
    swga.mkdirp(output_dir)

    kmers = []
    for k in xrange(min_size, max_size + 1):
        fg = primers.count_kmers(k, fg_genome_fp, output_dir, threshold)
        bg = primers.count_kmers(k, bg_genome_fp, output_dir, threshold)

        if exclude_fp:
            assert os.path.isfile(exclude_fp)
            ex = primers.count_kmers(k, exclude_fp, output_dir,
                                     exclude_threshold)
        else:
            ex = {}

        # Keep kmers found in foreground, merging bg binding values, and
        # excluding those found in the excluded fasta
        
        def primer_dict(i, seq):
            fg_freq = fg[seq]
            bg_freq = bg.get(seq, 0)
            ratio = fg_freq / float(bg_freq) if bg_freq > 0 else 0
            tm = Tm(seq)
            return {'seq': seq, 'fg_freq': fg_freq, 'bg_freq': bg_freq,
                    'ratio': ratio, 'tm': tm}
        
        kmers = (primer_dict(i, seq)
                 for i, seq in enumerate(fg.keys())
                 if seq not in ex.viewkeys())

        nkmers = len(fg.viewkeys() - ex.viewkeys())
        
        swga.message("Writing %d-mers into db..." % k)
        for kmer in progress.bar(kmers, expected_size=nkmers):
            Primer.create(**kmer)

    swga.message("Counted kmers in range %d-%d" % (min_size, max_size))
    
