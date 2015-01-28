'''
Provides methods for retrieving SWGA resources
'''
import yaml
from pkg_resources import (resource_exists, resource_filename, 
                           resource_stream)

def _get_resource_file(fp):
    resource_path = resource_filename("swga", fp)
    try:
        assert resource_exists("swga", fp)
        return resource_path
    except AssertionError:
        raise ValueError("The requested resource (%s) does not exist." 
                         % resource_path)


def get_dsk():
    dsk_bin = _get_resource_file("bin/dsk")
    parser_bin = _get_resource_file("bin/parse_results")
    return dsk_bin, parser_bin


def get_setfinder():
    sf_bin = _get_resource_file("bin/set_finder")
    return sf_bin


def get_swga_opts():
    with resource_stream("swga", "data/options.yaml") as opts_stream:
        opts = yaml.load(opts_stream)
        assert opts is not {}
        return opts


