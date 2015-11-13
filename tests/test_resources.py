import os
import swga.utils


def test_get_dsk():
    dsk_path = swga.utils.dsk
    assert dsk_path is not None and os.path.isfile(dsk_path)


def test_get_setfinder():
    sf_path = swga.utils.set_finder
    assert sf_path is not None and os.path.isfile(sf_path)
