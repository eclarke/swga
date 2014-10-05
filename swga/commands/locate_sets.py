import argparse
import swga
from pyfasta import Fasta
from swga import genome
import sys
import os

def main(argv, cfg_file, quiet):
    """
    Takes a list of primer sequences and returns the locations of each primer
    in the genome as a separate BED file for use in genome viewers.
    """
    parser = swga.basic_cmd_parser(description=main.__doc__,
                                   cmd_name="locate-sets",
                                   cfg_file=cfg_file)

    parser.add_argument('-i', '--input',
                        default=sys.stdin, metavar='FILE',
                        type=argparse.FileType('r'),
                        help="""File with primers sequences, one per line (default: stdin)""")

    parser.add_argument('-f', '--fasta',
                        required=True,
                        help="Path to fasta file (*not* flattened)")

    parser.add_argument('-o', '--output_folder',
                        required=True,
                        help="Path to folder to store BED files (one per primer)")
    
    args = parser.parse_args(argv)

    primers = [_.strip('\n') for _ in args.input.readlines()]
    
    fasta = Fasta(args.fasta)
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    for primer in primers:
        with open(os.path.join(args.output_folder, primer+'.bed'), 'w') as bedfile:
            bedfile.write("track name={}\n".format(primer))
            for record in fasta:
                seq = fasta[record][::]
                record_name = record.split('|')[0].strip()
                for location in genome.find_locations(primer, seq):
                    bedfile.write("{rec} {start} {end}\n".format(rec=record_name, 
                                                               start=location,
                                                               end=location+len(primer)))
                    
                
