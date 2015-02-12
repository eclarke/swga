import ConfigParser
import swga.utils.options
import pytest
import shlex
import sys
from StringIO import StringIO
from collections import OrderedDict


@pytest.fixture
def opts():
    opts = OrderedDict({
        'command': {
            'META': {
                'desc': 'test description',
            },
            'opt': {
                'desc': 'option description',
                'default': 10,
                'type': 'int',
            },
            'nocfg_opt': {
                'desc': 'option description',
                'incfg': False,
                'default': 'abc',
                'type': 'str'
            },
            'output': {
                'incfg': False,
                'default': 'stdout'
            }
        }
    })
    return opts
        

def test_cfg_from_opts(opts):
    cfg = StringIO(swga.utils.options.cfg_from_opts(opts))
    parser = ConfigParser.SafeConfigParser()
    parser.readfp(cfg)
    assert 'command' in parser.sections()
    assert 'opt' in dict(parser.items('command'))
    assert 'nocfg_opt' not in dict(parser.items('command'))


def test_argparser_from_opts(opts):
    parser = swga.utils.options.argparser_from_opts(opts, 'command')
    args = parser.parse_args(shlex.split("--nocfg_opt 20 --opt 5"))
    assert args.opt == 5
    assert args.nocfg_opt == "20"
    assert args.output == sys.stdout
