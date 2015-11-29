import os
import time
import signal
import subprocess

import utils


def find(**kwargs):
    """Find sets using the set_finder.

    :param min_bg_bind_dist: min avg bg binding distance between primers
    :param bg_length: bg genome length
    :param min_size: smallest allowable set size
    :param max_size: largest allowable set size
    :param graph_fp: compatibility graph file
    :param workers: if <= 1, searches graph in order of primer binding. If >1,
    searches graph from x randomly-chosen points in parallel.
    """
    workers = kwargs.get('workers', 1)
    if workers <= 1:
        return _find_sets(**{k:v for k,v in kwargs.items() if k !='workers'})
    else:
        return _mp_find_sets(**kwargs)


def _find_sets(min_bg_bind_dist, bg_length, min_size, max_size, graph_fp,
               vertex_ordering="weighted-coloring"):
    assert vertex_ordering in ["weighted-coloring", "random"]
    find_set_cmd = [
        utils.set_finder(), '-q', '-q',
        '--bg_freq', min_bg_bind_dist,
        '--bg_len', bg_length,
        '--min', min_size,
        '--max', max_size,
        '--unweighted',
        '--all',
        '--reorder', vertex_ordering,
        graph_fp
    ]
    find_set_cmd = " ".join([str(_) for _ in find_set_cmd])

    # We call the set_finder command as a line-buffered subprocess that
    # passes its output back to this process.
    # The function then yields each line as a generator; when close() is
    # called, it terminates the set_finder subprocess.
    process = subprocess.Popen(
        find_set_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        preexec_fn=os.setsid,
        bufsize=1)
    try:
        for line in iter(process.stdout.readline, b''):
            (yield line)
    finally:
        time.sleep(0.1)
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)


def _mp_find_sets(workers, **kwargs):
    setfinder_procs = [_find_sets(vertex_ordering="random", **kwargs)
                       for _ in range(workers)]
    try:
        while True:
            for i, setfinder in enumerate(setfinder_procs):
                (yield setfinder.next())
    finally:
        for i, setfinder in enumerate(setfinder_procs):
            setfinder.close()