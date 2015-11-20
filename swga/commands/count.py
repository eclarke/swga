# -*- coding: utf-8 -*-
import os
from collections import defaultdict

from swga import (error, message)
import swga.database as database
import swga.kmers
from swga.primers import Primer
from swga.commands._command import Command
from swga.utils import mkdirp
from peewee import OperationalError

INF = float('inf')
output_dir = ".swga_tmp"


class Count(Command):

    def run(self):
        if self.input:
            kmers = swga.kmers.parse_kmer_file(self.input)
            self.count_specific_kmers(kmers)
        else:
            self.count_kmers()

    def count_specific_kmers(self, kmers):
        try:
            # Skip primers that already exist and warn users
            existing = [p.seq for p in Primer.select().where(Primer.seq << kmers)]
            for p in existing:
                message("{} already exists in db, skipping...".format(p))
            kmers = [p for p in kmers if p not in existing]
        except OperationalError:
            # If this fails due to an OperationalError, it probably means the
            # database tables haven't been created yet. 
            error(
                "It doesn't appear that the workspace has been initialized: "
                "run `swga init' first.")
        mkdirp(output_dir)

        # Group the kmers by length to avoid repeatedly counting kmers of the
        # same size
        kmers_by_length = defaultdict(list)
        for kmer in kmers:
            kmers_by_length[len(kmer)].append(kmer)

        for k, mers in kmers_by_length.items():
            fg = swga.kmers.count_kmers(k, self.fg_genome_fp, output_dir, 1)
            bg = swga.kmers.count_kmers(k, self.bg_genome_fp, output_dir, 1)
            primers = []
            for mer in mers:
                try:
                    primers.append(primer_dict(mer, fg, bg, 0, INF, INF))
                except KeyError:
                    message(
                        "{} does not exist in foreground genome, skipping..."
                        .format(mer))

            # Omitting any primers that were returned empty
            # primers = filter(lambda p: p == {}, primers)
            chunk_size = 199
            message(
                "Writing {n} {k}-mers into db in blocks of {cs}..."
                .format(n=len(primers), k=k, cs=chunk_size))
            database.add_primers(primers, chunk_size, add_revcomp=False)

    def count_kmers(self):
        mkdirp(output_dir)

        kmers = []
        for k in xrange(self.min_size, self.max_size + 1):
            fg = swga.kmers.count_kmers(k, self.fg_genome_fp, output_dir)
            bg = swga.kmers.count_kmers(k, self.bg_genome_fp, output_dir)

            if self.exclude_fp:
                assert os.path.isfile(self.exclude_fp)
                ex = swga.kmers.count_kmers(
                    k, self.exclude_fp, output_dir, self.exclude_threshold)
            else:
                ex = {}

            # Keep kmers found in foreground, merging bg binding values, and
            # excluding those found in the excluded fasta

            kmers = [
                primer_dict(seq, fg, bg, self.min_fg_bind, self.max_bg_bind,
                            self.max_dimer_bp)
                for seq in fg.viewkeys() if seq not in ex.viewkeys()
            ]

            kmers = filter(lambda x: x != {}, kmers)

            nkmers = len(kmers)

            chunk_size = 199
            message(
                "Writing {n} {k}-mers into db in blocks of {cs}..."
                .format(n=nkmers * 2, k=k, cs=chunk_size))
            database.add_primers(kmers, chunk_size, add_revcomp=True)

        message("Counted kmers in range %d-%d" % (self.min_size, self.max_size))


def primer_dict(seq, fg, bg, min_fg_bind, max_bg_bind, max_dimer_bp):
    fg_freq = fg[seq]
    bg_freq = bg.get(seq, 0)
    ratio = fg_freq / float(bg_freq) if bg_freq > 0 else float('inf')

    # Setting to -1 disables the background binding frequency check
    max_bg_bind = INF if max_bg_bind < 0 else max_bg_bind

    if ((fg_freq >= min_fg_bind) and
            (bg_freq <= max_bg_bind) and
            # Homodimer check
            (swga.kmers.max_sequential_nt(seq, seq) <= max_dimer_bp)):
        return {
            'seq': seq,
            'fg_freq': fg_freq,
            'bg_freq': bg_freq,
            'ratio': ratio}
    else:
        return {}
