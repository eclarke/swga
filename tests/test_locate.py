import pytest
import swga.primers
import swga.locate
from swga.database import Primer


@pytest.fixture
def kmer():
    return "AAGG"


def test_locate_genome_positions(kmer, fastafile):
    locations = swga.locate.binding_sites(kmer, fastafile)
    assert len(locations['record1']) == 4
    assert locations['record2'] == []


def test_locate_chromosome_ends(fastafile):
    ends = swga.locate.chromosome_ends(fastafile)
    assert ends['record1'] == [0, 15]
    assert ends['record2'] == [16, 23]
    assert ends['record3'] == [24, 31]


def test_linearize_binding_sites(kmer, initdb, fastafile):
    p = Primer.create(seq=kmer)
    p._update_locations(fastafile)
    chr_ends = swga.locate.chromosome_ends(fastafile)
    linear_bind_sites = swga.locate.linearize_binding_sites([p], chr_ends)
    # (number of sites + (2*number of chromosomes) - (any overlaps))
    assert len(linear_bind_sites) == 10
    for record, ends in chr_ends.iteritems():
        start, end = ends
        assert start in linear_bind_sites
        assert end in linear_bind_sites
        for site in p.locations[record]:
            assert site in linear_bind_sites


def test_revcomp():
    assert "ATGC" == swga.locate.revcomp("GCAT")
