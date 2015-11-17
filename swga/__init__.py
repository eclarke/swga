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
    text = quote(msg, quote="!> ", nl=False)
    click.secho(text, err=True, fg='red')


def message(msg, newline=True):
    click.echo(msg, err=True, nl=newline)


def quote(text, width=72, quote="", nl=True):
    if not text:
        return ""
    out = ""
    for line in text.split("\n"):
        sublines = textwrap.wrap(line, width=width, replace_whitespace=False)
        sublines = [quote + l for l in sublines]# if l.strip()]
        out += "\n".join(sublines) + "\n"
    return out


class SWGAError(Exception):
    pass
