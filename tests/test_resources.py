import os
from collections import OrderedDict
import swga.utils.resources

def test_get_dsk():
    dsk_path = swga.utils.resources.get_dsk()
    assert dsk_path is not None and os.path.isfile(dsk_path)

def test_get_setfinder():
    sf_path = swga.utils.resources.get_setfinder()
    assert sf_path is not None and os.path.isfile(sf_path)

def test_get_swga_opts():
    assert isinstance(swga.utils.resources.get_swga_opts(), OrderedDict)
