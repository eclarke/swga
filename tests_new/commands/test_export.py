import os
from subprocess import check_call, check_output

import pytest


@pytest.mark.usefixtures('temporary_directory')
class TestExportCommands:

    def test_setup(self, simple_fastas):
        """Initialize the workspace for later export tests."""
        files = simple_fastas._asdict()
        print simple_fastas.fg
        print simple_fastas.bg
        check_call("swga init -f {fg} -b {bg} -e {ex}".format(**files), shell=True)
        check_call("swga count --input {primers}".format(**files), shell=True)
        check_call("swga activate {primers}".format(**files), shell=True)
        check_call("swga score --force --input {primers}".format(**files), shell=True)

    def test_lorenz(self):
        lorenz = check_output(
            "swga export lorenz --id -1 --no_header",
            shell=True
        )
        print lorenz
        # Retrieve the first line of the output and parse the results
        line = lorenz.split("\n")[0].split("\t")
        assert line[0] == "-1"

        # The results should all be floats and normalized to 0-1
        lz = [float(i) for i in line[1].split(",")]
        assert max(lz) == 1

    def test_bedgraph(self):
        """Bedgraphs should reflect binding sites, not nucleotides bound."""
        check_call(
            "swga export bedgraph --window_size=8 --step_size=8 --id -1",
            shell=True)
        bedgraph = "simple_fg_genome_export/set_-1/set_-1.bedgraph"
        assert os.path.isfile(bedgraph)
        with open(bedgraph, 'r') as bed:
            bed.readline()
            line2 = bed.readline()
            # If 16 was the recorded value for this interval, it's counting
            # nucleotides instead of primer hits
            assert line2 == "simple_genome 4 4 2\n"

    def test_bedfile(self):
        check_call(
            "swga export bedfile --id=-1", shell=True
        )
        bedfile = "simple_fg_genome_export/set_-1/AAAATTTT.bed"
        assert os.path.isfile(bedfile)
        with open(bedfile, 'r') as bed:
            assert bed.readline() == "track name=AAAATTTT\n"
            assert bed.readline() == "simple_genome 0 8\n"

