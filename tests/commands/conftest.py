import os
import shutil
import shlex
import tempfile
import pytest
from swga.commands import init


@pytest.fixture(scope="class")
def class_scoped_tmpdir(request):
    tmpdir = tempfile.mkdtemp()
    request.addfinalizer(lambda: shutil.rmtree(tmpdir))
    return tmpdir


@pytest.fixture(scope='class')
def workspace(simple_fastas, class_scoped_tmpdir, request):
    os.chdir(class_scoped_tmpdir)

    # Initialize a workspace in the tempdir
    argv = "-f {fg} -b {bg} -e {ex} --force".format(**simple_fastas._asdict())
    assert os.listdir(class_scoped_tmpdir) == []
    init.main(shlex.split(argv))
