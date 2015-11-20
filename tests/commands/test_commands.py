
import shlex
import pytest


#@pytest.mark.usefixtures('workspace')
@pytest.mark.xfail
class TestNormalWorkflow:

    def test_count(self):
        from swga.commands import Count
        argv = "--min_size 6 --max_size 6"
        c = Count(shlex.split(argv))
        c.run()

    def test_filter(self):
        from swga.commands import Filter
        argv = '--max_gini=1 --min_tm 0'
        c = Filter(shlex.split(argv))
        c.run()

    def test_find_sets(self):
        from swga.commands import FindSets
        c = FindSets([])
        c.run()

    def test_score(self):
        pass

    def test_export(self):
        pass

    def test_summary(self):
        pass

    def test_activate(self):
        pass
