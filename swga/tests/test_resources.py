import os
import pytest
import swga.utils


def test_get_dsk():
    """Must be able to find dsk."""
    dsk_path = swga.utils.dsk()
    assert dsk_path is not None and os.path.isfile(dsk_path)


def test_get_setfinder():
    """Must be able to find set_finder."""
    sf_path = swga.utils.set_finder()
    assert sf_path is not None and os.path.isfile(sf_path)


def test_specfiles():
    """Must be able to find all command specfiles."""
    from swga.commands import _commands
    for cmdname in _commands:
        swga.utils.specfile(cmdname)
