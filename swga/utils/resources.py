'''
Provides methods for retrieving SWGA resources
'''
import os
import sys
import yaml
import swga
from collections import OrderedDict
from pkg_resources import (
    resource_exists, resource_filename, resource_stream
)


def _get_resource_file(rs):
    _rs = os.path.join('bin', rs)
    res_path = os.path.join(sys.prefix, _rs)
    if not os.path.isfile(res_path) and resource_exists("swga", _rs):
        res_path = resource_filename("swga", _rs)
    else:
        swga.error("Could not find `{}': try re-installing swga.".format(rs))
    return res_path


def get_dsk():
    """Returns a path to the dsk executable."""
    return _get_resource_file('dsk')

def get_setfinder():
    """Returns a path to the set_finder executable."""
    return _get_resource_file('set_finder')


def get_swga_opts():
    with resource_stream("swga", "data/options.yaml") as opts_stream:
        opts = yaml.load(opts_stream)
        if opts == {}:
            swga.error("Empty options.yaml file: try reinstalling swga.")
        return opts

# This ensures that the options we get back from the yaml file are in order by
# replacing the default constructor with an OrderedDict
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG


def _dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())


def _dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_constructor(_mapping_tag, _dict_constructor)
yaml.add_representer(OrderedDict, _dict_representer)


