#import swga
import os
from pkg_resources import resource_stream
from ConfigParser import SafeConfigParser
import argutils
from argutils import (read, export)
from swga.clint.textui import indent, colored


class Command:
    args = None
    name = ""
    cfg_file = None

    def __init__(self, name, description=None, cfg_file=None):
        
        fp = os.path.join('commands', 'specfiles', name+'.yaml')
        specfile = resource_stream("swga", fp)
        opts = read.from_yaml(get_cmd_specfile(name))

        self.name = name
        self.parser = export.to_argparser(name, opts)
        config = SafeConfigParser()
        if config.read(cfg_file):
            parser = argutils.set_parser_defaults(name, opts)

    def parse_args(self, argv, quiet=False):
        args, unknown = self.parser.parse_known_args(argv)
        self.unknown = unknown
        self.args = vars(args)
        self.kwargs_as_args(**self.args)
        if not quiet:
            self.pprint_args()
        if len(self.args) > 0 and all(v is None for k, v in self.args.items()):
            swga.warn(
                "[swga %s]: All parameters are missing- ",
                "this may indicate a corrupt or missing parameters file."
                % self.name)

    def kwargs_as_args(self, **kwargs):
        '''
        Another interface to the command besides command-line arguments:
        specify them as arguments to this function.
        '''
        for k, v in kwargs.items():
            setattr(self, k, v)

    def pprint_args(self):
        if not self.args:
            return
        else:
            swga.message(colored.green("Command: " + self.name))
            with indent(2, quote='-'):
                for arg, val in self.args.iteritems():
                    if val is None or val == "":
                        val = colored.red("None")
                    swga.message(colored.blue("{}: {}".format(arg, val)))
