from __future__ import absolute_import
import click
import pkg_resources
import sys
import textwrap

DEFAULT_DB_FNAME = ".swga_data"
DEFAULT_CFG_FNAME = "parameters.cfg"

__version__ = pkg_resources.require("swga")[0].version


def error(msg, exception=True, wrap=True):
    '''Prints an error message to stderr and exits.'''
    if exception:
        raise SWGAError(msg)
    else:
        msg = quote(msg, width=75) if wrap else msg
        click.secho(msg, fg='red', err=True)
        sys.exit(1)


def warn(msg):
    '''Prints a warning message to stderr.'''
    text = quote(msg, quote="!> ")
    click.secho(text, err=True, fg='red')


def message(msg, newline=True):
    click.echo(msg, err=True, nl=newline)


def quote(text, width=72, quote="", nl=True):
    if not text:
        return ""
    lines = textwrap.wrap(text, width=width)
    lines = [quote + line.strip() for line in lines if line.strip()]
    out = "\n".join(lines)
    if nl:
        out += "\n"
    return out


class SWGAError(Exception):
    pass
