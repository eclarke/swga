from operator import attrgetter
from swga.primers import read_primer_file, write_primer_file
from swga.commands import Command


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)
    filter_primers(**cmd.args)


def filter_primers(input, 
                   output,
                   primers, 
                   max_bg_binding, 
                   num_primers):
    primers = read_primer_file(input)
    # sort by bg binding count
    primers = sorted(primers, key=attrgetter("bg_freq"))
    # remove primers that bind too many times to bg
    primers = [p for p in primers if p.bg_freq <= max_bg_binding]
    # keep only the top <num_primers>
    primers = primers[0:num_primers]
    # sort by fg/bg ratio
    primers = sorted(primers, key=attrgetter("ratio"), reverse=True)
    write_primer_file(primers, output)


def filter_melting_temp():
    pass


