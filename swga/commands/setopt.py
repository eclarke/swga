import click
import swga
from ConfigParser import SafeConfigParser
from swga.commands import Command


def main(argv, cfg_file):
    cmd = Command('setopt', cfg_file=cfg_file)
    cmd.parse_args(argv, quiet=True)
    setopt(cfg_file, **cmd.args)


def setopt(cfg_file, command, opt, value):
    parser = SafeConfigParser()
    cfg = parser.read(cfg_file)
    if not cfg:
        swga.error("No config file found at %s" % cfg_file)

    if not parser.has_section(command):
        swga.error("Specified command section [%s] does not exist in %s" % (command, cfg_file))

    old_val = None
    if parser.has_option(command, opt):
        old_val = parser.get(command, opt)
    click.confirm("Change [{command}]:{opt} from {old_val} --> {value}?".format(**locals()), abort=True)
    parser.set(command, opt, value)
    with open(cfg_file, 'wb') as out:
        parser.write(out)
