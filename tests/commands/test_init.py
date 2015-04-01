import os
from click.testing import CliRunner

from swga.commands.init import (
    main,
    DEFAULT_PARAMETERS_FNAME)

def __init_testing_partial(input_list):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(main, input = "\n".join(input_list))
        assert result.exit_code == 0
        assert os.path.isfile(DEFAULT_PARAMETERS_FNAME)

        
def test_init_noprompts(fastafile):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            main, ['-f', fastafile, '-b', fastafile, '-e', fastafile])
        print result.output
        assert result.exit_code == 0
        assert os.path.isfile(DEFAULT_PARAMETERS_FNAME)

        
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
