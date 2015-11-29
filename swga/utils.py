import os
import errno
from pkg_resources import (
    resource_exists, resource_filename, resource_stream
)
import functools

import swga
from click import progressbar

__all__ = [
    "mkdirp",
    "chunk_iterator",
    "specfile",
    "dsk",
    "set_finder"
]




def mkdirp(path):
    """Create a directory unless it already exists, using EAFP methods."""
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def chunk_iterator(itr, fn, n=100, show_progress=True, label=None):
    """Break an iterable into chunks and apply a function to each chunk.

    :param itr: the iterable to be chunked
    :param fn: the function to be applied to each chunks
    :param n: the size of each chunk
    :param show_progress: show a progress bar
    :param label: the label to show on the progress bar
    """
    length = len(itr)
    if length == 0:
        return
    label = "" if label is None else label
    if length / n <= 1:
        show_progress = False
    chunked_itr = chunks(itr, n)
    if show_progress:
        with progressbar(chunked_itr, max(length / n, 1), label) as bar:
            for chunk in bar:
                fn(chunk)
    else:
        for chunk in chunked_itr:
            fn(chunk)


def _get_resource_file(rs):
    import sys
    _rs = os.path.join(sys.prefix, 'bin', rs)
    # If it's not in sys.prefix/bin/, try sys.exec_prefix?
    if not os.path.isfile(_rs):
        _rs = os.path.join(sys.exec_prefix, 'bin', rs)
    # If it still doesn't work, check the package data
    if not os.path.isfile(_rs):
        if resource_exists('swga', os.path.join('bin', rs)):
            _rs = resource_filename('swga', os.path.join('bin', rs))
        else:
            swga.error("Could not find `{}': try re-installing swga.".format(rs))
    return os.path.abspath(_rs)


def specfile(name):
    """Return the path to the specfile for a given command name."""
    fp = os.path.join('commands', 'specfiles', name + '.yaml')
    return resource_stream("swga", fp)


def fmtkv(k, v):
    """Pretty-print a key-value pair."""
    set_stat_line = "{key:.<16}: {val: >10s}"
    num_fmt_str = "{:,G}"
    if not isinstance(v, basestring):
        print k, v
        v = num_fmt_str.format(v)
    return set_stat_line.format(key=k, val=v)


def dsk():
    return _get_resource_file('dsk')


def set_finder():
    return _get_resource_file('set_finder')
