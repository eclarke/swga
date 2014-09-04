import argparse
import subprocess

def main(argv, cfg_file, quiet):
    '''
    Flatten a FASTA file into one line, stripping entry headers and
    keeping only the > entry markers.
    '''
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='flatten',
                                   cfg_file=cfg_file)
    
    parser.add_argument('-i', '--input', metavar="FILE"
                        required=True, 
                        help="Input FASTA file")

    parser.add_argument('-o', '--output', metavar="FILE",
                        type=argparse.FileType('w', 0),
                        help="Output filename (default:
                        input_name.flattened)")
    
    args = parser.parse_args(argv)
    if not os.path.isfile(args.input):
        sys.stderr.write("Error: input file specified does not "
                         "exist!\n")
        exit(1)
                         
    if not args.output:
        args.output = args.input+'.flattened'
    
    cmdstr = "sed 's/>.*/>/' {input} | tr -d '\n' > {output}"
    subprocess.check_call(cmdstr.format(**vars(args)))
    if not quiet:
        sys.stderr.write("{input} --> {output}".format(**vars(args)))
    
        

