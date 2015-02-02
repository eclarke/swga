from swga.commands2 import Command
import swga.primers as primers


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_kmers(**cmd.args)


def filter_kmers(fg_genome_fp,
                 bg_genome_fp,
                 kmer_dir,
                 min_avg_fg_bind,
                 max_avg_bg_bind,
                 min_tm,
                 max_tm,
                 num_primers):
    
                 
