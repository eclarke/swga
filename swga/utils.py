import os
import sys
import errno
import swga
from click import progressbar
from pkg_resources import (
    resource_exists, resource_filename
)


def mkdirp(path):
    '''Simulates 'mkdir -p': creates a directory unless it already exists'''
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def chunk_iterator(itr, fn, n=100, show_progress=True, label=None):
    '''Breaks an iterable into chunks and applies a function to each chunk.
    Arguments:
    - itr the iterable to be chunked
    - fn the function to be applied to each chunks
    - n the size of each chunk
    - show_progress show a progress bar
    - label the label to show on the progress bar
    '''
    length = len(itr)
    if length == 0:
        return
    label = "" if label is None else label
    if length/n <= 1:
        show_progress = False
        swga.message(label)
    chunked_itr = chunks(itr, n)
    if show_progress:
        with progressbar(chunked_itr, max(length/n, 1), label) as bar:
            for chunk in bar:
                fn(chunk)
    else:
        for chunk in chunked_itr:
            fn(chunk)


def _get_resource_file(rs):
    _rs = os.path.join('bin', rs)
    #res_path = os.path.join(sys.prefix, _rs)
    if resource_exists("swga", _rs):
        res_path = resource_filename("swga", _rs)
        return res_path
    else:
        swga.error("Could not find `{}': try re-installing swga.".format(rs))


dsk = _get_resource_file('dsk')
set_finder = _get_resource_file('set_finder')
