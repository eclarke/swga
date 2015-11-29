import pytest
from swga import workspace
from swga.workspace import Set, Primer
from swga.primers import Primers

class TestWorkspace(object):
    def test_workspace_connection(self):
        """Workspace should be initialized."""
        import sys
        with workspace.connection(':memory:') as ws:
            assert not ws.is_closed()
            assert ws.database == ':memory:'

    def test_create_tables(self):
        with workspace.connection(':memory:') as ws:
            ws.create_tables()
            tables = ws.get_tables()
            assert len(tables) == len(workspace._tables)

    def test_workspace_metadata(self):
        """Metadata getting and setting should use property attributes."""
        with workspace.connection(':memory:') as ws:
            ws.create_tables()
            ws.metadata = {
                'version': '1.0.0',
                'fg_file': 'somewhere'
            }
            m = ws.metadata
            assert m.version == '1.0.0'
            assert m.fg_file == 'somewhere'
            assert m.bg_file is None


@pytest.mark.usefixtures('ws')
class TestSets(object):

    def test_set_add(self, primers):
        # Test expected behavior: set added if it doesn't exist
        Set.add(0, primers, score=1)
        assert Set.get(Set._id == 0)
        # Duplicate sets should not be added to the database
        _, created = Set.add(1, primers, score=1)
        assert not created

    def test_bad_set_add(self):
        with pytest.raises(ValueError):
            Set.add(0, None, score=1)
        with pytest.raises(ValueError):
            invalid_primers = Primer.select().where(Primer.seq == "XX")
            Set.add(0, primers=invalid_primers, score=1)

@pytest.mark.usefixtures('ws')
class TestPrimers(object):

    def test_add_primers(self):
        '''Must add the reverse complement of a primer if requested.'''
        primers = [{'seq': "AAAA"}]
        Primers.add(primers, add_revcomp=True)
        assert Primer.select().where(Primer.seq == "TTTT").count() == 1
