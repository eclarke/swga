import pytest

import swga.database
from swga.database import Primer, Set, PrimerSet


@pytest.fixture
def seqs():
    return ['A' * (i + 1) for i in range(20)]

@pytest.fixture(scope="function")
def initdb(request):
    swga.database.db.init(":memory:")
    swga.database.db.create_tables([Primer, Set, PrimerSet])
    def fin():
        print ("Closing database")
        swga.database.db.drop_tables([Primer, Set, PrimerSet])
        swga.database.db.close()
    request.addfinalizer(fin)

@pytest.fixture
def tprimers(seqs):
    return [Primer.create(seq=seq) for seq in seqs]

@pytest.fixture
def tset():
    return Set.create(_id=1, score=100) 


