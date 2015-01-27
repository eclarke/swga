from swga import parse_config
from swga.utils.options import load_swga_opts, argparser_from_opts

import count, filter, flatten, init, locate_mers, mkgraph, score, sets

class Command:

    def __init__(self, name, description=None, cfg_file=None):
        opts = load_swga_opts()
        self.parser = argparser_from_opts(opts, name)

        defaults = parse_config(cfg_file, name) if cfg_file else {}
        self.parser.set_defaults(**defaults)
        

    def parse_args(self, argv):
        self.args = vars(self.parser.parse_args(argv))
        self.kwargs_as_args(**self.args)
        
        
    def kwargs_as_args(self, **kwargs):
        '''
        Another interface to the command besides command-line arguments:
        specify them as arguments to this function. The values will replace the
        command-line ones if conflicting, and it is not necessary to call
        self.parse_args beforehand.
        '''
        for k, v in kwargs.items():
            setattr(self, k, v)



    def run(self, **kwargs):
        pass
