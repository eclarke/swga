from __future__ import absolute_import
from .clint.textui import puts, colored, STDERR, indent, max_width
import pkg_resources
import sys
import textwrap

__version__ = pkg_resources.require("swga")[0].version


def error(msg, exception=True):
    '''Prints an error message to stderr and exits.'''
    if exception:
        raise SWGAError(msg)
    else:
        _errprint(msg)
        sys.exit(1)


def warn(msg):
    '''Prints a warning message to stderr.'''
    with indent(3, quote=colored.red("!> ")):
        _errprint(msg)


def message(msg, newline=True):
    puts(msg, stream=STDERR, newline=newline)


def _errprint(text):
    text = colored.red(max_width(textwrap.dedent(text), 75))
    puts(text, stream=STDERR)


class SWGAError(Exception):
    pass
