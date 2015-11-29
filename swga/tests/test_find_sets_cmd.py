import pytest
import shlex
import swga.score
import swga.sets
from swga.workspace import Primer
from swga.graph import build_graph

@pytest.fixture
def primers():
    primers = [
        Primer.create(seq="ATGC", fg_freq=1, bg_freq=2, ratio=1.0, active=True),
        Primer.create(seq="GGCC", fg_freq=1, bg_freq=3, ratio=0.5, active=True),
        Primer.create(seq="CCTA", fg_freq=2, bg_freq=0, ratio=float('inf'), active=True)
    ]
    return primers
    

def test_build_graph(ws, primers, tmpdir):
    outfile = tmpdir.join('graph')
    build_graph(max_hetdimer_bind=3, outfile=str(outfile))
    graph = outfile.read()
    print graph
    assert graph == """p sp 3 3
n 1 1
n 2 2
n 3 3
e 1 2
e 1 3
e 2 3
"""


def test_find_sets(tmpdir):
    from swga.commands import FindSets
    graph = ('p sp 6 7\n'
             'n 1 1\n'
             'n 2 2\n'
             'n 3 1\n'
             'n 4 4\n'
             'n 5 5\n'
             'n 6 6\n'
             'e 1 2\n'
             'e 1 3\n'
             'e 2 3\n'
             'e 2 4\n'
             'e 4 5\n'
             'e 4 6\n'
             'e 5 6\n')
    fp = tmpdir.join("testgraph")
    fp.write(graph)
    sets = swga.sets.find(
        min_bg_bind_dist=2,
        min_size=3,
        max_size=3,
        bg_length=10,
        graph_fp=fp)
    output = list(sets)
    assert len(output) == 1
    ids, bg_dist_mean = swga.score.read_set_finder_line(output[0])
    assert ids == [1, 2, 3]
    assert bg_dist_mean == 2.5


