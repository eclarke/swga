import os
import errno
import argparse
import subprocess
import multiprocessing
import swga

def main(argv, cfg_file, quiet):
    """
    Count k-mers in a foreground genome, then count those k-mers in
    a larger background genome.
    """
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='count',
                                   cfg_file=cfg_file)

    parser.add_argument('-f', '--fg_genome', 
                        metavar="FASTA",
                        required=True,
                        help='''FASTA file containing the foreground
                        or target genome (default: %(default)s)''')

    parser.add_argument('-b', '--bg_genome', 
                        metavar="FASTA",
                        required=True,
                        help='''FASTA file containing the background
                        genome (default: %(default)s)''')

    parser.add_argument('-m', '--min', 
                        type=int, 
                        help="""Min primer size (default:
                        %(default)s)""")

    parser.add_argument('-M', '--max', 
                        type=int, 
                        help="""Max primer size (default:
                        %(default)s)""")

    parser.add_argument('-o', '--output_dir', 
                        help="""Directory to store kmer FASTA
                        files.""")

    args = parser.parse_args(argv)

    if not os.path.isfile(args.fg_genome):
        parser.print_usage()
        parser.exit(1, ("Error: specified foreground genome file %s "
                        "does not exist" % args.fg_genome))

    if args.bg_genome is None or not os.path.isfile(args.bg_genome):
        parser.print_usage()
        parser.exit(1, ("Error: specified background genome file %s "
                        "does not exist" % args.bg_genome))

    if args.min > args.max:
        parser.print_usage()
        parser.exit(1, ("Error: Max primer size must be larger than min primer "
                        "size.\n"))

    count_mers(args)


def count_mers(args):
    fg_genome_name = args.fg_genome.split(os.sep).pop()
    bg_genome_name = args.bg_genome.split(os.sep).pop()
    fg_kmer_dir = os.path.join(args.output_dir, fg_genome_name+'.kmers')
    bg_kmer_dir = os.path.join(args.output_dir, bg_genome_name+'.kmers')
    # Try to create the directory, unless it already exists
    for d in (args.output_dir, fg_kmer_dir, bg_kmer_dir):
        try:
            os.makedirs(d)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    for i in xrange(args.min, args.max+1):
        fg_output_name = os.path.join(fg_kmer_dir, '%i-mers.fa' % i)
        bg_output_name = os.path.join(bg_kmer_dir, '%i-mers.fa' % i)
        fg_cmd_str = ("{jellyfish} --mer-len={i} --size=100M --output={output} "
            "--threads={threads} {genome}").format(
                jellyfish=args.jellyfish,
                i=i,
                threads=args.nthreads,
                output=fg_output_name,
                genome=args.fg_genome)
        bg_cmd_str = ("{jellyfish} --mer-len={i} --size=100M --output={output} "
            "--threads={threads} --if={fgp} {genome}").format(
                jellyfish=args.jellyfish,
                i=i,
                threads=args.nthreads,
                fgp=fg_output_name,
                output=bg_output_name,
                genome=args.bg_genome)
        print fg_cmd_str
        subprocess.check_call(fg_cmd_str, shell=True)
        print bg_cmd_str
        subprocess.check_call(bg_cmd_str, shell=True)



if __name__ == '__main__':
    main(None, None)
