import swga
import argparse
from swga.utils.options import load_swga_opts, add_args_from_opts

class Command:

    def __init__(self, name, description=None, cfg_file=None):
        self.parser = argparse.ArgumentParser(
            prog = 'swga ' + name,
            description = description)

        opts = load_swga_opts()
        add_args_from_opts(self.parser, name, opts)
        if cfg_file:
            defaults, _ = swga.parse_config(cfg_file, name)
        else:
            defaults = {}
        self.parser.set_defaults(**defaults)

        
    def parse_args(self, argv):
        self.args = self.parser.parse_args(argv)
    

    def run(self, **kwargs):
        pass
