import pytest
import swga.database as database
from swga.database import Primer, Set, Primer_Set
from playhouse.test_utils import test_database

class TestPrimersSets:
    
    @pytest.fixture
    def seqs(self):
        return ['A' * (i+1) for i in range(10)]

    @pytest.fixture
    def db(self):
        db = database.init_db(':memory:')
        return db

    def create_primers(self, seqs):
        for seq in seqs:
            Primer.create(seq=seq)

    def test_primers_sets(self, seqs, db):
        with test_database(db, (Primer, Set, Primer_Set), fail_silently=False):
            self.create_primers(seqs)
            p2 = Primer.create(seq="TTTT")
            seq_subset = seqs[0:5]
            primers = list(Primer.select().where(Primer.seq << seq_subset).execute())
            newset = Set.create(_id=1, score=100)
            for primer in primers:
                Primer_Set.create(seq=primer.seq, set=newset)
            assert newset.primers.count() == 5
            p_in_s = (Set
                      .select()
                      .join(Primer_Set)
                      .join(Primer)
                      .where(Set._id == newset._id))
            assert p_in_s.count() == 5
            for primer in database.get_primers_for_set(newset._id):
                assert primer in primers
                assert p2 not in primers
            
