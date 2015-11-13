import os
import click
import argutils
import argutils.read
import argutils.export
from pkg_resources import resource_stream
from ConfigParser import SafeConfigParser
from swga import (
    warn, message, DEFAULT_DB_FNAME, DEFAULT_CFG_FNAME,
    __version__
)
import swga.database as database


_commands = [
    'activate',
    'count',
    'export',
    'filter',
    'find_sets',
    'score',
    'summary',
]


def specfile(name):
    fp = os.path.join('commands', 'specfiles', name + '.yaml')
    return resource_stream("swga", fp)


def create_config_file():
    cfg_file_str = ""
    for cmd in _commands:
        spec = specfile(cmd)
        opts = argutils.read.from_yaml(spec)
        if opts:
            if opts['_meta'].get("_exclude"):
                continue
            print argutils.export.to_config(cmd, opts)
            cfg_file_str += argutils.export.to_config(cmd, opts) + "\n"
    return cfg_file_str


class Command:
    args = None
    name = ""
    cfg_file = None

    def __init__(self, name,
                 cfg_file=DEFAULT_CFG_FNAME,
                 db_name=DEFAULT_DB_FNAME):

        spec = specfile(name)
        opts = argutils.read.from_yaml(spec)

        self.name = name
        self.parser = argutils.export.to_argparser(name, opts)
        config = SafeConfigParser()
        if config.read(cfg_file):
            self.parser = argutils.set_parser_defaults(self.parser, config)

        database.init_db(db_name)
        database.check_version()
        meta = database.Metadata.get()

    def parse_args(self, argv, quiet=False):
        args, unknown = self.parser.parse_known_args(argv)
        self.unknown = unknown
        self.args = vars(args)
        self.kwargs_as_args(**self.args)
        if not quiet:
            self.pprint_args()
        if len(self.args) > 0 and all(v is None for k, v in self.args.items()):
            warn(
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
            message(click.style("swga v{}: {}".format(__version__, self.name), fg='green'))
            for arg, val in self.args.iteritems():
                if val is None or val == "":
                    val = click.style("None", fg='red')
                message(click.style("  - {}: {}".format(arg, val), fg='blue'))
