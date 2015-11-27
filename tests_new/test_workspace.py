from swga import workspace


class TestWorkspace(object):

    def test_workspace_connection(self):
        """Workspace should be initialized."""
        import sys
        with workspace.connection(':memory:') as ws:
            sys.stderr.write("Checking if closed...")
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
            assert m.bg_file == ''
