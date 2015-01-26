from operator import attrgetter
from swga.primers import read_primer_file, write_primer_file
from swga.commands2 import Command

class Filter(Command):
    
    def run(self):
        primers = read_primer_file(self.args['input'])
        primers = filter_primers(primers, 
                                 self.args['max_bg_binding'], 
                                 self.args['num_primers'])
        write_primer_file(primers, self.args['output'])


def filter_primers(primers, max_bg_binding, num_primers):
    # sort by bg binding count
    primers = sorted(primers, key=attrgetter("bg_freq"))
    # remove primers that bind too many times to bg
    primers = [p for p in primers if p.bg_freq <= max_bg_binding]
    # keep only the top <num_primers>
    primers = primers[0:num_primers]
    # sort by fg/bg ratio
    primers = sorted(primers, key=attrgetter("ratio"), reverse=True)
    return primers
