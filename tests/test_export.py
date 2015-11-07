import os
from subprocess import check_call

def test_setup(isolated_filesystem, simple_fastas):
    """Initialize the workspace for later export tests."""
    with isolated_filesystem:
        files = simple_fastas._asdict()
        print simple_fastas.fg
        print simple_fastas.bg
        check_call("swga init -f {fg} -b {bg} -e {ex}".format(**files), shell=True)
        check_call("swga count --input {primers}".format(**files), shell=True)
        check_call("swga activate {primers}".format(**files), shell=True)
        check_call("swga score --force --input {primers}".format(**files), shell=True)

def test_bedgraph(isolated_filesystem):
    """Bedgraphs should reflect binding sites, not nucleotides bound."""
    with isolated_filesystem:
        check_call(
            "swga export bedgraph --window_size=8 --step_size=8 --id -1", 
            shell=True)
        bedgraph = "simple_fg_genome_export/set_-1/set_-1.bedgraph"
        assert os.path.isfile(bedgraph)
        print open(bedgraph).read()
        with open(bedgraph, 'r') as bed:
            bed.readline()
            line2 = bed.readline()
            # If 16 was the recorded value for this interval, it's counting 
            # nucleotides instead of primer hits
            assert line2 == "simple_genome 4 4 2\n"


def test_bedfile(isolated_filesystem):
    with isolated_filesystem:
        check_call(
            "swga export bedfile --id=-1", shell=True
        )
        bedfile = "simple_fg_genome_export/set_-1/AAAATTTT.bed"
        assert os.path.isfile(bedfile)
        with open(bedfile, 'r') as bed:
            assert bed.readline() == "track name=AAAATTTT\n"
            assert bed.readline() == "simple_genome 0 8\n"

