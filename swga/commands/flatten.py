import subprocess
import swga
import sys
import os

def main(argv, cfg_file, quiet):
    '''
    Flatten a FASTA file into one line, stripping entry headers and
    keeping only the > entry markers.
    '''
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name='flatten',
                                   cfg_file=cfg_file)
    
    parser.add_argument('-i', '--input', metavar="FILE",
                        required=True, 
                        help="Input FASTA file")

    parser.add_argument('-o', '--output', metavar="FILE",
                        help="Output filename (default: \
                        input_name.flattened)")

    parser.add_argument('-f', '--force', action='store_true',
                        help="Overwrite output file if it already exists") 
    
    args = parser.parse_args(argv)
    if not os.path.isfile(args.input):
        sys.stderr.write("Error: input file specified does not "
                         "exist.\n")
        exit(1)
                         
    if not args.output:
        args.output = args.input+'.flattened'

    if os.path.isfile(args.output):
        if args.force:
            os.remove(args.output)
        else:
            sys.stderr.write("Error: output file exists.\n")
            exit(1)

    cmdstr = "sed 's/>.*/>/' {input} | tr -d '\\n' > {output}"
    cmdstr = cmdstr.format(**vars(args))
    subprocess.check_call(cmdstr.format(**vars(args)), shell=True)
    if not quiet:
        sys.stderr.write("{input} --> {output}\n".format(**vars(args)))
    
        

