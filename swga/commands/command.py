from swga.utils.options import argparser_from_opts
from swga.utils.resources import get_swga_opts
from swga import parse_config, message, warn
from swga.clint.textui import indent, colored

class Command:
    args = None
    name = ""
    cfg_file = None
    def __init__(self, name, description=None, cfg_file=None):
        opts = get_swga_opts()
        self.name = name
        self.parser = argparser_from_opts(opts, name)

        defaults = parse_config(cfg_file, name) if cfg_file else {}
        self.parser.set_defaults(**defaults)
        

    def parse_args(self, argv, quiet=False):
        args, unknown = self.parser.parse_known_args(argv)
        self.unknown = unknown
        self.args = vars(args)
        self.kwargs_as_args(**self.args)
        if not quiet:
            self.pprint_args()
        if len(self.args) > 0 and all(v is None for k, v in self.args.items()):
            warn("[swga %s]: All parameters are missing- this may indicate a corrupt or "
                 "missing parameters file." % self.name)
        
    def kwargs_as_args(self, **kwargs):
        '''
        Another interface to the command besides command-line arguments:
        specify them as arguments to this function. The values will replace the
        command-line ones if conflicting, and it is not necessary to call
        self.parse_args beforehand.
        '''
        for k, v in kwargs.items():
            setattr(self, k, v)

    def pprint_args(self):
        if not self.args: 
            return
        else:
            message(colored.green("Command: " + self.name))
            with indent(2, quote='-'):
                for arg, val in self.args.iteritems():
                    if val is None or val == "":
                        val = colored.red("None")
                    message(colored.blue("{}: {}".format(arg, val)))
