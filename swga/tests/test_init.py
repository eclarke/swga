import os
import __builtin__
import pytest
import click
from click.testing import CliRunner

from swga.commands import init
from swga import DEFAULT_CFG_FNAME

@pytest.fixture
def runner(monkeypatch):
    def run(input_list):
        monkeypatch.setattr(__builtin__, 'raw_input', '\n'.join(input_list))
        init.main([])
    return run


def test_init_noprompts(fastafile):
    with CliRunner().isolated_filesystem():
        init.main(['-f', fastafile, '-b', fastafile, '-e', fastafile])
        assert os.path.isfile(DEFAULT_CFG_FNAME)


def test_init_prompts(runner, fastafile):
    with CliRunner().isolated_filesystem():
        runner([
            fastafile,  # foreground
            fastafile,  # background
            "Y",        # yes, specify exclude_fp
            fastafile   # exclude_fp
        ])
        assert os.path.isfile(DEFAULT_CFG_FNAME)


def test_init_noexclude(runner, fastafile):
    with CliRunner().isolated_filesystem():
        runner([
            fastafile,
            fastafile,
            "N"  # no exclude_fp
        ])
        assert os.path.isfile(DEFAULT_CFG_FNAME)