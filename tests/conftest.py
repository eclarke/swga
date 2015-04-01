import os
import shutil
import tempfile
import contextlib
import pytest
import swga.database
from swga.database import Primer, Set, PrimerSet


@pytest.fixture
def seqs():
    return ['A' * (i + 1) for i in range(20)]

@pytest.fixture(scope="function")
def initdb(request):
    swga.database.db.init(":memory:")
    swga.database.db.create_tables([Primer, Set, PrimerSet])
    def fin():
        print ("Closing database")
        swga.database.db.drop_tables([Primer, Set, PrimerSet])
        swga.database.db.close()
    request.addfinalizer(fin)

@pytest.fixture
def tprimers(seqs):
    return [Primer.create(seq=seq) for seq in seqs]

@pytest.fixture
def tset():
    return Set.create(_id=1, score=100) 

@pytest.fixture
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

@pytest.fixture
def fg_fasta(request, testdata_fp):
    '''Foreground genome in integration tests'''
    return _fastafile(request, testdata_fp, "bburgdorferi_wgs.fa")

@pytest.fixture
def bg_fasta(request, testdata_fp):
    '''Background genome in integration tests'''
    return _fastafile(request, testdata_fp, "ecoli_wgs.fa")

@pytest.fixture
def ex_fasta(request, testdata_fp):
    '''Sequences to exclude in integration tests'''
    return _fastafile(request, testdata_fp, "bburgdorferi_plasmid_lp28-4.fa")


@pytest.fixture(scope="module")
def module_scoped_tmpdir(request):
    tmpdir = tempfile.mkdtemp()
    def fin():
        shutil.rmtree(tmpdir)
    request.addfinalizer(fin)
    return tmpdir


@pytest.fixture(scope="function")
def isolated_filesystem(module_scoped_tmpdir):

    @contextlib.contextmanager
    def tmpdir_context_manager(tmpdir):
        start_path = os.getcwd()
        os.chdir(tmpdir)        
        yield
        os.chdir(start_path)

    return tmpdir_context_manager(module_scoped_tmpdir)
