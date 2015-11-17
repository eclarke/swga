import os
import shlex
import pytest

# From http://stackoverflow.com/questions/12411431/pytest-how-to-skip-the-rest-of-tests-in-the-class-if-one-has-failed
# Enables the incremental testing of commands, and xfails the rest if one fails
# (since often the others will fail anyway for unhelpful reasons)
def pytest_runtest_makereport(item, call):
    if "incremental" in item.keywords:
        if call.excinfo is not None:
            parent = item.parent
            parent._previousfailed = item


def pytest_runtest_setup(item):
    previousfailed = getattr(item.parent, "_previousfailed", None)
    if previousfailed is not None:
        pytest.xfail("previous test failed (%s)" % previousfailed.name)


@pytest.mark.incremental
class TestCommands:

    def test_init(self, isolated_filesystem, simple_fastas):
        from swga.commands import init
        from swga import DEFAULT_DB_FNAME, DEFAULT_CFG_FNAME
        with isolated_filesystem:
            argv = "-f {fg} -b {bg} -e {ex}".format(**simple_fastas._asdict())
            init.main(shlex.split(argv))
            assert os.path.isfile(DEFAULT_CFG_FNAME)
            assert os.path.isfile(DEFAULT_DB_FNAME)

    def test_count(self, isolated_filesystem):
        from swga.commands import Count
        with isolated_filesystem:
            argv = "--min_size 6 --max_size 6"
            c = Count(shlex.split(argv))
            c.run()

    def test_filter(self, isolated_filesystem):
        pass

    def test_find_sets(self, isolated_filesystem):
        pass

    def test_score(self, isolated_filesystem):
        pass

    def test_export(self, isolated_filesystem):
        pass

    def test_summary(self, isolated_filesystem):
        pass

    def test_activate(self, isolated_filesystem):
        pass
