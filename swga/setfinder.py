from swga import message
import swga.utils.resources as resources
import subprocess
import os
import time
import signal


def find_sets(
        min_bg_bind_dist,
        min_size,
        max_size,
        bg_genome_len,
        graph_fp,
        vertex_ordering="weighted-coloring"):
    assert vertex_ordering in ["weighted-coloring", "random"]
    set_finder = resources.get_setfinder()

    find_set_cmd = [set_finder, '-q', '-q', '-B', min_bg_bind_dist,
                    '-L', bg_genome_len, '-m', min_size, '-M',
                    max_size, '-a', '-u', '-r', vertex_ordering,
                    graph_fp]
    find_set_cmd = " ".join([str(_) for _ in find_set_cmd])

    # We call the set_finder command as a line-buffered subprocess that passes
    # its output back to this process. The function then yields each line as a
    # generator; when close() is called, it terminates the set_finder
    # subprocess.
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


def mp_find_sets(nprocesses=4, **kwargs):
    '''
    Runs multiple independent copies of the set finder using randomized vertex
    coloring. This means each process explores a different solution space
    (though none are necessarily optimal) in hopes that good sets are found
    faster.
    '''
    setfinders = [find_sets(vertex_ordering='random', **kwargs)
                  for _ in xrange(nprocesses)]
    try:
        while True:
            for i, setfinder in enumerate(setfinders):
                (yield setfinder.next())
    finally:
        for i, setfinder in enumerate(setfinders):
            setfinder.close()
