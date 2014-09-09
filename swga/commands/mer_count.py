import swga
import subprocess

def main(argv, cfg_file, quiet):
    """
    Automates process of finding kmer counts for a range of k's. Uses genometools.
    """
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='count',
                                   cfg_file=cfg_file)

    parser.add_argument("-f", "--fg_genome",
                        metavar="fasta",
                        required=True,
                        help="Foreground or target genome (FASTA format)")

    parser.add_argument("-b", "--bg_genome",
                        metavar="fasta",
                        required=True,
                        help="Background genome (FASTA format)")

    parser.add_argument("-m", "--min", type=int,
                        metavar="min",
                        help="Minimum kmer size (default: %(default)s)")

    parser.add_argument("-M", "--max", type=int,
                        metavar="max",
                        help="Maximum kmer size (default: %(default)s)")

    parser.add_argument("-n", "--min_num", type=int,
                        metavar="num",
                        help="""Minimum # of times kmer must appear
                        (default: %(default)s)""")


    suffix_opts = parser.add_mutually_exclusive_group()

    suffix_opts.add_argument("--esa",
                             metavar="FILE",
                             help="""Name of existing suffix array (if
                             not specified, creates one)""")

    args = parser.parse_args(argv)
    # create the suffix array if it does not already exist for the specified range
    if not args.esa:
        args.esa = args.fg_genome+'.esa.%i-%i'%(args.min, args.max)
        make_esa(args)
    # count mers in fg and bg genome for given range
    count_mers(args, quiet)


def make_esa(args):
    '''
    Creates an enhanced suffix array (esa) for the foreground
    genome.
    '''
    mk_esa_cmd = ("gt suffixerator -dna -pl -tis -suf -lcp -v "
                  "-parts 4 -db {fg_genome} "
                  "-indexname {esa}").format(**vars(args))
    print ">> "+mk_esa_cmd
    subprocess.check_call(mk_esa_cmd, shell=True)
    

def count_mers(args, quiet):
    '''
    Iteratively counts k-mers in the foreground and background genome
    along a range of values of k.
    '''
    mkindex_cmdstr = ("gt tallymer mkindex -scan -counts -pl "
                      "-mersize {i} "
                      "-minocc {min_num} "
                      "-esa {esa} "
                      "-indexname {idxname}")
    search_cmdstr = ("gt tallymer search -output counts sequence "
                     "-tyr {idxname} -q {query_genome} > {output}")

    for i in xrange(args.min, args.max+1):
        idxname = args.fg_genome+'.idx.%i' % i
        fg_output = args.fg_genome+'.%imers' % i
        bg_output = args.fg_genome+'.%imers' % i
        mkindex_cmd = mkindex_cmdstr.format(
            i=i,
            idxname=idxname,
            min_num=args.min_num,
            esa=args.esa)
        fg_search_cmd = search_cmdstr.format(
            idxname=idxname,
            query_genome=args.fg_genome,
            output=fg_output)
        bg_search_cmd = search_cmdstr.format(
            idxname=idxname,
            query_genome=args.bg_genome,
            output=bg_output)
        # build fg index
        print ">> "+mkindex_cmd
        subprocess.check_call(mkindex_cmd, shell=True)
        # find fg bind #s
        print ">> "+fg_search_cmd
        subprocess.check_call(fg_search_cmd, shell=True)
        # find bg bind #s
        print ">> "+bg_search_cmd
        subprocess.check_call(bg_search_cmd, shell=True)
