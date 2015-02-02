from swga.commands2 import Command
import swga.primers as primers
from swga.primers_pw import Primer, db
import os
import swga


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_kmers(**cmd.args)


def filter_kmers(kmer_dir,
                 min_avg_fg_bind,
                 max_avg_bg_bind,
                 min_tm,
                 max_tm,
                 num_primers):
    if not os.path.isfile(os.path.join(kmer_dir, 'primers.db')):
        swga.swga_error("Missing primers database: re-run `swga count`")

