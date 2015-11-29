import os
import pytest
import swga
from swga import DEFAULT_DB_FNAME as db_fname
import swga.workspace
from swga.workspace import Primer, Set


def pytest_addoption(parser):
    parser.addoption('--integration', action='store_true', help='Run integration tests')


@pytest.fixture(scope='function')
def ws(tmpdir):
    db_fn = str(tmpdir.join(db_fname))
    with swga.workspace.connection(db_fn) as _ws:
        _ws.create_tables()


@pytest.fixture
def seqs():
    return ['A' * (i + 1) for i in range(20)]


@pytest.fixture
def primers(ws, seqs):
    return [Primer.create(seq=seq) for seq in seqs]


@pytest.fixture
def tset(ws):
    return Set.create(_id=1, score=100)


@pytest.fixture(scope='session')
def testdata_fp():
    return os.path.join(os.path.dirname(__file__), "data")


def _fastafile(request, testdata_fp, fname):
    '''Returns a given FASTA file from tests/data, cleaning up afterward'''
    fa = os.path.join(testdata_fp, fname)
    fai = fa + ".fai"

    def fin():
        if os.path.exists(fai):
            os.remove(fai)
    request.addfinalizer(fin)
    return fa


@pytest.fixture
def fastafile(request, testdata_fp):
    '''Fake fasta file for location tests'''
    return _fastafile(request, testdata_fp, "test.fasta")


@pytest.fixture(scope='module')
def simple_fastas(request, testdata_fp):
    from collections import namedtuple
    Files = namedtuple('Files', ['fg', 'bg', 'ex', 'primers'])
    fg = _fastafile(request, testdata_fp, "simple_fg_genome.fa")
    bg = _fastafile(request, testdata_fp, "simple_bg_genome.fa")
    ex = _fastafile(request, testdata_fp, "simple_ex_seqs.fa")
    primers = _fastafile(request, testdata_fp, "simple_primers.txt")
    return Files(fg, bg, ex, primers)


@pytest.fixture(scope='class')
def temporary_directory(request):
    import tempfile, shutil
    tmpdir = tempfile.mkdtemp()
    curdir = os.curdir
    os.chdir(tmpdir)

    def fin():
        os.chdir(curdir)
        shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def metadata():
    swga.meta(
        db_name='',
        version=swga.__version__,
        fg_file='',
        bg_file='',
        ex_file='',
        fg_length=0,
        bg_length=0
    )