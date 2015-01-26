from swga import parse_config
from swga.utils.options import load_swga_opts, mk_argparser_from_opts

class Command:

    def __init__(self, name, description=None, cfg_file=None):

        opts = load_swga_opts()
        self.parser = mk_argparser_from_opts(opts, name)

        defaults , _ = parse_config(cfg_file, name) if cfg_file else {}, None
        self.parser.set_defaults(**defaults)
        
        self.args = None
        

    def parse_args(self, argv):
        self.args = vars(self.parser.parse_args(argv))
        
        
    def kwargs_as_args(self, **kwargs):
        '''
        Another interface to the command besides command-line arguments:
        specify them as arguments to this function. The values will replace the
        command-line ones if conflicting, and it is not necessary to call
        self.parse_args beforehand.
        '''

        if not self.args:
            self.args = kwargs
            return

        for k, v in kwargs:
            if k in self.args:
                self.args[k] = v


    def run(self, **kwargs):
        pass
