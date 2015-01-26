import pytest
from swga.core import parse_config


class TestParseConfig:

    @pytest.fixture()
    def config_file(self, tmpdir):
        cfg_file = tmpdir.join("test.cfg")
        cfg_file.write("""[section1]
val1 = 1
val2 = 
val3 = 3
""")
        return cfg_file

    def test_empty_config(self, config_file):
        """Parser should accept missing sections."""
        defaults, cfg_parser = parse_config(str(config_file), "missing_section")
        assert "missing_section" not in cfg_parser.sections()
        assert defaults == {}

    def test_missing_vals(self, config_file):
        """Parser should accept missing values."""
        defaults, cfg_parser = parse_config(str(config_file), "section1")
        assert defaults == {"val1":'1', "val2":'', "val3":'3'}
