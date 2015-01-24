import sys
from operator import attrgetter
from swga.primers import read_primer_file
from swga.commands2 import Command


class Filter(Command):
    
    def run(self, **kwargs):
        primerline = "{seq}\t{fg_freq}\t{bg_freq}\t{ratio}\n"
        primers = self.filter_primers(**kwargs)
        for primer in primers:
            self.args.output.write(primerline.format(**primer._asdict()))


    def filter_primers(self, input_fp, max_bg_binding, num_primers):
        primers = read_primer_file(input_fp, False)
        # sort by bg binding count
        primers = sorted(primers, key=attrgetter("bg_freq"))
        # remove primers that bind too many times to bg
        primers = [p for p in primers if p.bg_freq <= max_bg_binding]
        # keep only the top <num_primers>
        primers = primers[0:num_primers]
        # sort by fg/bg ratio
        primers = sorted(primers, key=attrgetter("ratio"), reverse=True)
        return primers


def test_filter():
    cmd = Filter('filter', 'description')
    cmd.parse_args(['--input', 'infile',
                    '--output', 'outfile',
                    '--max_bg_binding', 5000,
                    '--num_primers', 200])
    assert cmd.args == {'input': 'infile',
                        'output': 'outfile',
                        'max_bg_binding': 5000,
                        'num_primers': 200}
