# -*- coding: utf-8 -*-
import os
import swga
import swga.primers
import swga.database as database
from collections import defaultdict
from swga.core import chunk_iterator
from swga.primers import Primer
from swga.commands import Command
import click

output_dir = ".swga_tmp"

def main(argv, cfg_file):
    cmd = Command('count', cfg_file=cfg_file)
    cmd.parse_args(argv)
    database.init_db(cmd.primer_db, create_if_missing=True)
    if cmd.input:
        kmers = swga.primers.parse_kmer_file(cmd.input)
        count_specific_kmers(kmers, cmd.fg_genome_fp, cmd.bg_genome_fp)
    else:
        count_kmers(**cmd.args)


def count_specific_kmers(kmers, fg_genome_fp, bg_genome_fp):
    # Remove primers that already exist and warn users
    existing = [p.seq for p in Primer.select().where(Primer.seq << kmers)]
    for p in existing:
        swga.message("{} already exists in db, skipping...".format(p))
    kmers = filter(lambda p: p not in existing, kmers)

    # Group the kmers by length to avoid repeatedly counting kmers of the same size
    kmers_by_length = defaultdict(list)
    for kmer in kmers:
        kmers_by_length[len(kmer)].append(kmer)

    for k, mers in kmers_by_length.items():
        fg = swga.primers.count_kmers(k, fg_genome_fp, output_dir, 1)            
        bg = swga.primers.count_kmers(k, bg_genome_fp, output_dir, 1)
        primers = [primer_dict(mer, fg, bg, 0, float('inf')) for mer in mers]
        chunk_size = 199
        swga.message("Writing {n} {k}-mers into db in blocks of {cs}..."
                     .format(n=len(primers), k=k, cs=chunk_size))
        chunk_iterator(primers,
                       fn=lambda c: Primer.insert_many(c).execute(),
                       n=chunk_size, label="Updating database...")
    
    
def count_kmers(fg_genome_fp,
                bg_genome_fp, 
                min_size, 
                max_size, 
                min_fg_bind,
                max_bg_bind,
                primer_db,
                exclude_fp,
                exclude_threshold, **kwargs):
    assert os.path.isfile(fg_genome_fp)
    assert os.path.isfile(bg_genome_fp)

    if os.path.isfile(primer_db):
        swga.warn("Existing database found at %s" % os.path.abspath(primer_db))
        swga.warn("Re-counting primers will reset the entire database!")
        click.confirm("Are you sure you want to proceed?", abort=True)

    database.init_db(primer_db, create_if_missing=True)
    database.create_tables()
    

    swga.mkdirp(output_dir)

    kmers = []
    for k in xrange(min_size, max_size + 1):
        fg = swga.primers.count_kmers(k, fg_genome_fp, output_dir, 1)
        bg = swga.primers.count_kmers(k, bg_genome_fp, output_dir, 1)

        if exclude_fp:
            assert os.path.isfile(exclude_fp)
            ex = swga.primers.count_kmers(k, exclude_fp, output_dir,
                                     exclude_threshold)
        else:
            ex = {}

        # Keep kmers found in foreground, merging bg binding values, and
        # excluding those found in the excluded fasta
        
        kmers = [primer_dict(seq, fg, bg, min_fg_bind, max_bg_bind)
                 for seq in fg.viewkeys()
                 if seq not in ex.viewkeys()]

        kmers = filter(lambda x: x != {}, kmers)
        
        nkmers = len(kmers)

        chunk_size = 199
        swga.message("Writing {n} {k}-mers into db in blocks of {cs}..."
                     .format(n=nkmers, k=k, cs=chunk_size))
        chunk_iterator(kmers,
                       fn=lambda c: Primer.insert_many(c).execute(),
                       n=chunk_size, label="Blocks written: ")

    swga.message("Counted kmers in range %d-%d" % (min_size, max_size))
    

def primer_dict(seq, fg, bg, min_fg_bind, max_bg_bind):
    fg_freq = fg[seq]
    bg_freq = bg.get(seq, 0) 
    ratio = fg_freq / float(bg_freq) if bg_freq > 0 else 0
    if fg_freq < min_fg_bind or bg_freq > max_bg_bind:
        return {}
    return {'seq': seq, 'fg_freq': fg_freq, 'bg_freq': bg_freq,
            'ratio': ratio}
