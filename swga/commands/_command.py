from ConfigParser import SafeConfigParser
import abc
import click
import argutils
import swga.database as database
from swga import (
    DEFAULT_CFG_FNAME,
    DEFAULT_DB_FNAME,
    __version__,
    warn,
    message
)
import swga.utils as utils


class Command(object):
    args = None
    name = ""
    cfg_file = None

    def __init__(
            self,
            name,
            cfg_file=DEFAULT_CFG_FNAME,
            db_name=DEFAULT_DB_FNAME):

        spec = utils.specfile(name)
        opts = argutils.read.from_yaml(spec)

        self.name = name
        self.parser = argutils.export.to_argparser(name, opts)
        config = SafeConfigParser()
        if config.read(cfg_file) and config.has_section(name):
            self.parser = argutils.set_parser_defaults(self.parser, config)

        database.init_db(db_name)
        database.check_version()
        meta = database.Metadata.get()
        self.primer_db = db_name
        self.fg_genome_fp = meta.fg_file
        self.bg_genome_fp = meta.bg_file
        self.exclude_fp = meta.ex_file
        self.fg_length = meta.fg_length
        self.bg_length = meta.bg_length

    @abc.abstractmethod
    def run(self):
        '''Run the command.'''
        return

    def parse_args(self, argv, quiet=False):
        args, unknown = self.parser.parse_known_args(argv)
        self.unknown = unknown
        self.args = vars(args)
        self._kwargs_as_args(**self.args)
        if not quiet:
            self.pprint_args()
        if len(self.args) > 0 and all(v is None for k, v in self.args.items()):
            warn(
                "[swga %s]: All parameters are missing- ",
                "this may indicate a corrupt or missing parameters file."
                % self.name)

    def _kwargs_as_args(self, **kwargs):
        '''Sets command-line arguments as object attributes.'''
        for k, v in kwargs.items():
            setattr(self, k, v)

    def pprint_args(self):
        if not self.args:
            return
        else:
            message(
                click.style("swga {}, v{}".format(self.name, __version__),
                            fg='green'))
            for arg, val in self.args.iteritems():
                if val is None or val == "":
                    val = click.style("None", fg='red')
                message(click.style("  - {}: {}".format(arg, val), fg='blue'))
