import os
import json
import pytest


import swga.locate
from swga.database import Primer

@pytest.fixture()
def fastafile(request):
    fastafile = os.path.join(os.path.dirname(__file__), "data", "test.fasta")
    fastaidx = fastafile + '.fai'
    # remove the index after each test so that it's recreated
    def fin():
        os.remove(fastaidx)
    request.addfinalizer(fin)
    return fastafile

@pytest.fixture
def kmer():
    return "ATGC"

 
def test_locate_genome_positions(kmer, fastafile):
    locations = swga.locate.binding_sites(kmer, fastafile)
    assert len(locations['record1']) == 32
    assert locations['record2'] == []


def test_locate_chromosome_ends(fastafile):
    ends = swga.locate.chromosome_ends(fastafile)
    assert ends['record1'] == [0, 127]
    assert ends['record2'] == [128, 191]
    assert ends['record3'] == [192, 199]

    
def test_linearize_binding_sites(kmer, initdb, fastafile):
    p = Primer.create(seq=kmer)
    p.locations = json.dumps(swga.locate.binding_sites(p.seq, fastafile))
    chr_ends = swga.locate.chromosome_ends(fastafile)
    linear_bind_sites = swga.locate.linearize_binding_sites([p], chr_ends)
    # (number of sites + (2*number of chromosomes) - (any overlaps))
    assert len(linear_bind_sites) == 38
    for record, ends in chr_ends.iteritems():
        start, end = ends
        assert start in linear_bind_sites
        assert end in linear_bind_sites
        for site in json.loads(p.locations)[record]:
            assert site in linear_bind_sites
        
