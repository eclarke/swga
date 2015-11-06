import os
import sys
import errno
import swga
from swga.clint.textui import progress


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
    if len(itr) == 0:
        return
    label = "" if label is None else label
    if len(itr)/n <= 1:
        show_progress = False
        swga.message(label)
    chunked_itr = chunks(itr, n)
    if show_progress:
        chunked = progress.bar(
            chunked_itr,
            label=label,
            expected_size=max(len(itr)/n, 1))
    else:
        chunked = chunked_itr
    for chunk in chunked:
        fn(chunk)
