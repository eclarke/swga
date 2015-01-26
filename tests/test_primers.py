import pytest
from swga.primers import Primer, read_primer_file, write_primer_file
class TestPrimers:

    @pytest.fixture()
    def primer_file(self, tmpdir):
        f = tmpdir.join("primers.tmp")
        f.write("""AAAA 1 2 3
# some malformed line here
GGGG 1 2 3
CCCC 2 3 4 4""")
        return f

    @pytest.fixture()
    def primer_list(self):
        pl = [Primer(1, "AAAA", fg_freq=1, bg_freq=2, ratio=3.0),
              Primer(3, "GGGG", fg_freq=1, bg_freq=2, ratio=3.0),
              Primer(4, "CCCC", fg_freq=2, bg_freq=3, ratio=4.0)]
        return pl

    def test_read_primer_file(self, primer_file, primer_list):
        primers = read_primer_file(primer_file.readlines())
        assert primers
        print primers
        assert all([primer in primer_list for primer in primers])

    def test_write_primer_file(self, primer_list, tmpdir):
        tmpout = tmpdir.join("test_primers.tmp")
        write_primer_file(primer_list, tmpout)
        primers = read_primer_file(tmpout.readlines())
        assert primers
        print primers
        for primer in primers:
            assert primer in primer_list

