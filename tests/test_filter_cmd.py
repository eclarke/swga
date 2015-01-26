import pytest
from StringIO import StringIO
from swga.commands2.filter import Filter
from swga.primers import Primer, read_primer_file

@pytest.fixture(scope="module")
def primer_file():
    f = StringIO("""Some header file here
ATGC 50 99 0.5
TCGA 100 200 0.4
CGTA\t100\t98\t0.2
GCCT;1;2;3
CTTA 4 99 0.4
""")
    return f

@pytest.fixture(scope="module")
def primer_list():
    pl = [Primer(2, "ATGC", fg_freq=50, bg_freq=99, ratio=0.5),
          Primer(4, "CGTA", fg_freq=100, bg_freq=98, ratio=0.2)]
    return pl


@pytest.fixture(scope="module")
def output_file():
    return StringIO()


@pytest.fixture(scope="module")
def filter_cmd(primer_file, output_file):
    filter_cmd = Filter('filter')
    filter_cmd.kwargs_as_args(input=primer_file,
                              output=output_file,
                              max_bg_binding=200,
                              num_primers=2)
    return filter_cmd


def test_filter_cmd(filter_cmd, output_file, primer_list):
    '''
    Filter should remove invalid primers and reorder the two valid remaining
    primers.
    '''
    filter_cmd.run()
    primers = read_primer_file(output_file.getvalue().split("\n"))
    assert primers
    for primer in primers:
        print primer
        assert primer in primer_list
