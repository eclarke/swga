# -*- coding: utf-8 -*-
import os
import swga
import swga.primers as primers
from swga.primers_pw import Primer, db
from swga.commands2 import Command


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

    swga.mkdirp(output_dir)

    db.init(os.path.join(output_dir, 'primers.db'))
    db.connect()
    db.drop_tables([Primer], safe=True)
    db.create_tables([Primer])

    kmers = []    
    for k in xrange(min_size, max_size + 1):
        fg = primers.count_kmers(k, fg_genome_fp, output_dir, threshold)
        bg = primers.count_kmers(k, bg_genome_fp, output_dir, threshold)
        
        kmers = [primers.Primer(id=i, seq=seq, fg_freq=fg[seq], 
                                bg_freq=bg.get(seq, 0))
                 for i, seq in enumerate(fg.keys())]
        
        kmers = ({'seq': seq, 'fg_freq': fg[seq], 'bg_freq':bg.get(seq, 0), 
                  'ratio': fg[seq]/float(bg.get(seq, 0))} 
                 for i, seq in enumerate(fg.keys()))
        swga.message("Writing %d-mers into db..." % k)
        for kmer in kmers:
            Primer.create(**kmer)

    swga.message("Counted kmers in range %d-%d" % (min_size, max_size))
    
