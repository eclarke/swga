"""
Integration tests for swga.
"""
import os
from subprocess import check_call, check_output
import pytest


@pytest.fixture
def fastafile(request, testdata_fp):
    def fasta(fname):
        fa = os.path.join(testdata_fp, fname)
        fai = fa + '.fai'

        def rm_fai():
            try:
                os.remove(fai)
            except OSError:
                pass

        request.addfinalizer(rm_fai)

        return fa
    return fasta


@pytest.fixture
def fg_fasta(fastafile):
    """Foreground genome in integration tests"""
    return fastafile("bburgdorferi_wgs.fa")


@pytest.fixture
def bg_fasta(fastafile):
    """Background genome in integration tests"""
    return fastafile("ecoli_wgs.fa")


@pytest.fixture
def ex_fasta(fastafile):
    """Sequences to exclude in integration tests"""
    return fastafile("bburgdorferi_plasmid_lp28-4.fa")


@pytest.mark.integration
@pytest.mark.usefixtures('temporary_directory')
class TestIntegration:

    def test_init(self, fg_fasta, bg_fasta, ex_fasta):
        command = "swga init -f {} -b {} -e {} --force".format(
            fg_fasta, bg_fasta, ex_fasta)
        check_call(command, shell=True)

    def test_count_repeat(self):
        """Repeat count commands should reset primers."""
        command = "swga count --min_size 5 --max_size 5 --min_tm 0 --max_tm 100"
        check_call(command, shell=True)
        command = "swga count --min_size 5 --max_size 5 --min_tm 0 --max_tm 100 --force"
        check_call(command, shell=True)

    def test_count(self):
        command = "swga count --min_size 7 --max_size 7 --max_bg_bind -1 --force"
        check_call(command, shell=True)

    def test_filter(self):
        command = "swga filter --max_bg_bind 1000000 --min_tm 0 --max_tm 100"
        check_call(command, shell=True)

    def test_find_sets(self):
        command = "swga find_sets"
        check_call(command, shell=True)

    def test_find_sets_randomized(self):
        command = "swga find_sets --workers=2 --force"
        check_call(command, shell=True)

    def test_export_sets(self):
        command = "swga export sets --limit 1 --order_by score"
        check_call(command, shell=True)

    def test_export_bedfile(self):
        command = "swga export bedfile --id 3"
        check_call(command, shell=True)

    def test_export_bedgraph(self):
        command = "swga export bedgraph --id 3"
        check_call(command, shell=True)

    def test_count_manual(self):
        with open('primers.txt', 'w') as out:
            out.write("ATGCAT\nATTTAT\n")
        command = "swga count --input primers.txt"
        check_call(command, shell=True)

    def test_activate(self):
        command = "swga activate primers.txt"
        check_call(command, shell=True)

    def test_summary(self):
        command = "swga summary"
        check_call(command, shell=True)


@pytest.mark.integration
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