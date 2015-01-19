import os
import subprocess
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

    swgahome = swga.get_swgahome()
    args.kmer_counter = os.path.join(swgahome, 'bin', 'kmer_locations')
    if not quiet:
        swga.print_cfg_file(parser.prog, cfg_file)
        swga.print_args(parser.prog, args)
    count_mers(args)


def count_mers(args):
    fg_genome_name = args.fg_genome.split(os.sep).pop()
    bg_genome_name = args.bg_genome.split(os.sep).pop()
    fg_kmer_dir = os.path.join(args.output_dir, fg_genome_name+'.fgkmers')
    bg_kmer_dir = os.path.join(args.output_dir, bg_genome_name+'.bgkmers')

    # Create the directories if they don't exist already
    for d in (args.output_dir, fg_kmer_dir, bg_kmer_dir):
        swga.mkdirp(d)

    for i in xrange(args.min, args.max+1):
        fg_output_name = os.path.join(fg_kmer_dir, '%i-mers.fa' % i)
        bg_output_name = os.path.join(bg_kmer_dir, '%i-mers.fa' % i)
        fg_cmd_str = "{kmer_counter} --kmer {i} --input {genome} --label > {output}".format(
            kmer_counter=args.kmer_counter,
            i=i,
            genome=args.fg_genome,
            output=fg_output_name)
        fg_mers_only = fg_output_name+".mersonly"
        bg_cmd_str = ("{kmer_counter} --kmer {i} --input {genome} --mer-file {fgmers} "
                      "--label > {output}").format(
                          kmer_counter=args.kmer_counter,
                          i=i,
                          genome=args.bg_genome,
                          fgmers=fg_mers_only,
                          output=bg_output_name)
        print fg_cmd_str
        subprocess.check_call(fg_cmd_str, shell=True)
        # copy all the mer seqs to a temporary file (since kmer_counter expects only seqs, not also counts)
        subprocess.check_call("awk '{print $1}' < %s > %s" % (fg_output_name, fg_mers_only), shell=True)
        print bg_cmd_str
        subprocess.check_call(bg_cmd_str, shell=True)
        os.remove(fg_mers_only)



if __name__ == '__main__':
    main(None, None)
