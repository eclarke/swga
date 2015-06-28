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


def progressbar(i, length):
    if i >= 1:
        i = i/(length*1.0)
    sys.stderr.write('\r[%-20s] %-3d%%' % ('='*int(round(i*20)), i*100))
    sys.stderr.flush()


def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def chunk_iterator(itr, fn, n=100, show_progress=True, label=None):
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
