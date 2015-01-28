import swga
import argparse
from swga.utils.options import get_swga_opts, add_args_from_opts

class Command:

    def __init__(self, name, description=None, cfg_file=None):
        self.parser = argparse.ArgumentParser(
            prog = 'swga ' + name,
            description = description)

        opts = get_swga_opts()
        add_args_from_opts(self.parser, name, opts)
        defaults, _ = swga.parse_config(cfg_file, name)
        self.parser.set_defaults(**defaults)

        
    def parse_args(self, argv):
        self.args = self.parser.parse_args(argv)
    

    def run(self, **kwargs):
        pass
