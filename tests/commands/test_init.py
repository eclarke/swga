import os
from click.testing import CliRunner

from swga.commands import init
from swga import DEFAULT_CFG_FNAME


def __init_testing_partial(input_list):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(init.main, input="\n".join(input_list))
        assert result.exit_code == 0
        assert os.path.isfile(DEFAULT_CFG_FNAME)


def test_init_noprompts(fastafile):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            init.main, ['-f', fastafile, '-b', fastafile, '-e', fastafile])
        print result.output
        assert result.exit_code == 0
        assert os.path.isfile(DEFAULT_CFG_FNAME)


def test_init_prompts(fastafile):
    __init_testing_partial([
        fastafile,  # foreground
        fastafile,  # background
        "Y",        # yes, specify exclude_fp
        fastafile   # exclude_fp
    ])


def test_init_noexclude(fastafile):
    __init_testing_partial([
        fastafile,
        fastafile,
        "N"  # no exclude_fp
    ])
