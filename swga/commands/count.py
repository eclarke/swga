# -*- coding: utf-8 -*-
import os
import swga
import swga.primers
import swga.database as database
from collections import defaultdict
from swga.primers import Primer, max_consecutive_binding
from swga.commands import Command
from peewee import OperationalError
import click

INF = float('inf')

output_dir = ".swga_tmp"

def main(argv, cfg_file):
    cmd = Command('count', cfg_file=cfg_file)
    cmd.parse_args(argv)
    database.init_db(cmd.primer_db, create_if_missing=True)
    if cmd.input:
        kmers = swga.primers.parse_kmer_file(cmd.input)
        count_specific_kmers(kmers, **cmd.args)
    else:
        count_kmers(**cmd.args)

        
def check_create_tables(primer_db):
    if os.path.isfile(primer_db):
        swga.warn("Existing database found at %s" % os.path.abspath(primer_db))
        swga.warn("This will reset the entire database!")
        click.confirm("Are you sure you want to proceed?", abort=True)
    database.create_tables()

        
def count_specific_kmers(
        kmers,
        fg_genome_fp,
        bg_genome_fp,
        primer_db,
        **kwargs):
    
    try:
        # Skip primers that already exist and warn users
        existing = [p.seq for p in Primer.select().where(Primer.seq << kmers)]
        for p in existing:
            swga.message("{} already exists in db, skipping...".format(p))
        kmers = filter(lambda p: p not in existing, kmers)
    except OperationalError:
        # If this fails due to an OperationalError, it probably means the
        # database tables haven't been created yet
        check_create_tables(primer_db)
        swga.mkdirp(output_dir)
        
    # Group the kmers by length to avoid repeatedly counting kmers of the same size
    kmers_by_length = defaultdict(list)
    for kmer in kmers:
        kmers_by_length[len(kmer)].append(kmer)

    for k, mers in kmers_by_length.items():
        fg = swga.primers.count_kmers(k, fg_genome_fp, output_dir, 1)            
        bg = swga.primers.count_kmers(k, bg_genome_fp, output_dir, 1)
        primers = []
        for mer in mers:
            try:
                primers.append(primer_dict(mer, fg, bg, 0, INF, INF))
            except KeyError:
                swga.message(
                    "{} does not exist in foreground genome, skipping..."
                    .format(mer)) 
        
        # Omitting any primers that were returned empty
        primers = filter(lambda p: p == {}, primers)
        chunk_size = 199
        swga.message(
            "Writing {n} {k}-mers into db in blocks of {cs}..."
            .format(n=len(primers), k=k, cs=chunk_size))
        database.add_primers(primers, chunk_size, add_revcomp=False)

        
def count_kmers(
        fg_genome_fp,
        bg_genome_fp, 
        min_size, 
        max_size, 
        min_fg_bind,
        max_bg_bind,
        max_dimer_bp,
        primer_db,
        exclude_fp,
        exclude_threshold,
        **kwargs):
    assert os.path.isfile(fg_genome_fp)
    assert os.path.isfile(bg_genome_fp)

    check_create_tables(primer_db)
    swga.mkdirp(output_dir)

    kmers = []
    for k in xrange(min_size, max_size + 1):
        fg = swga.primers.count_kmers(k, fg_genome_fp, output_dir)
        bg = swga.primers.count_kmers(k, bg_genome_fp, output_dir)

        if exclude_fp:
            assert os.path.isfile(exclude_fp)
            ex = swga.primers.count_kmers(k, exclude_fp, output_dir,
                                     exclude_threshold)
        else:
            ex = {}

        # Keep kmers found in foreground, merging bg binding values, and
        # excluding those found in the excluded fasta
        
        kmers = [primer_dict(seq, fg, bg, min_fg_bind, max_bg_bind, max_dimer_bp)
                 for seq in fg.viewkeys()
                 if seq not in ex.viewkeys()]

        kmers = filter(lambda x: x != {}, kmers)
        
        nkmers = len(kmers)

        chunk_size = 199
        swga.message("Writing {n} {k}-mers into db in blocks of {cs}..."
                     .format(n=nkmers*2, k=k, cs=chunk_size))
        database.add_primers(kmers, chunk_size, add_revcomp=True)

    swga.message("Counted kmers in range %d-%d" % (min_size, max_size))
    

def primer_dict(seq, fg, bg, min_fg_bind, max_bg_bind, max_dimer_bp):
    fg_freq = fg[seq]
    bg_freq = bg.get(seq, 0) 
    ratio = fg_freq / float(bg_freq) if bg_freq > 0 else float('inf')
    if ((fg_freq >= min_fg_bind) and
        (bg_freq <= max_bg_bind) and
        (max_consecutive_binding(seq, seq) <= max_dimer_bp)):
        return {
            'seq': seq,
            'fg_freq': fg_freq,
            'bg_freq': bg_freq,
            'ratio': ratio}
    else:
        return {}
