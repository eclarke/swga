import pytest
import swga
import swga.core
import swga.database as database
from swga.database import Primer, Set

class TestPrimersSets:
    
    def test_create_primers_sets(self, initdb, tprimers, tset):
        '''
        Create primers and sets, and find which primers belong to which sets.
        '''
        mers = Primer.select().limit(5)
        tset.primers.add(mers)
        assert mers.count() == 5
        for mer in mers:
            assert tset in mer.sets
            assert mer in tset.primers

    def test_destroy_primer_sets(self, initdb, tprimers, tset):
        '''
        The many-to-many relationship must be updated when a set is destroyed.
        '''
        mers = Primer.select().limit(5)
        tset.primers.add(mers)
        tset.delete_instance()
        assert mers.count() == 5
        for mer in mers:
            assert tset not in mer.sets

    def test_add_set_function(self, initdb, tprimers):
        '''Tests adding a set using add_set function.'''
        s = database.add_set(_id=2, primers=tprimers, score=100)
        assert s.primers.count() == len(tprimers)
        for primer in s.primers:
            assert primer in tprimers

    def test_bad_add_set_function(self, initdb, tprimers):
        '''Should raise errors if invalid primers supplied.'''
        with pytest.raises(swga.core.SWGAError):
            database.add_set(_id=2, primers=None, score=100)
        with pytest.raises(swga.core.SWGAError):
            invalid_primers = Primer.select().where(Primer.seq == "XX")
            database.add_set(_id=3, primers=invalid_primers, score=100)

    def test_update_in_chunks(self, initdb, tprimers, seqs):
        '''Must push all the updates successfully..'''
        for primer in tprimers:
            primer.fg_freq = 100
        database.update_in_chunks(tprimers)
        primers = Primer.select().where(Primer.seq << seqs)
        for primer in primers:
            assert primer.fg_freq == 100

    def test_create_tables(self):
        database.db.init(":memory:")
        database.create_tables()
        p = Primer.create(seq="ATGC")
        s = Set.create(_id=1, score=1)
        s.primers.add(p)
        database.db.close()
        
    def test_add_primers(self, initdb):
        '''Must add the reverse complement of a primer if requested.'''
        primers = [{'seq': "AAAA"}]
        database.add_primers(primers, add_revcomp=True)
        assert Primer.select().where(Primer.seq == "TTTT").count() == 1
    
