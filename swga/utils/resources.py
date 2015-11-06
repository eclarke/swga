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


def _get_resource_file(fp):
    resource_path = resource_filename("swga", fp)
    if not resource_exists("swga", fp):
        swga.error("Resource does not exist: {}".format(resource_path))
    return resource_path


def get_dsk():
    dsk = os.path.join(sys.prefix, 'bin', 'dsk')
    if not os.path.isfile(dsk):
        swga.error("Cannot find `dsk' in `{}'; try reinstalling swga."
                   .format(dsk))
    return dsk


def get_setfinder():
    setfinder = os.path.join(sys.prefix, 'bin', 'set_finder')
    if not os.path.isfile(setfinder):
        swga.error("Cannot find `set_finder' in `{}'; try reinstalling swga."
                   .format(setfinder))
    return setfinder


def get_swga_opts():
    with resource_stream("swga", "data/options.yaml") as opts_stream:
        opts = yaml.load(opts_stream)
        if opts == {}:
            swga.swga_error("Empty options.yaml file: try reinstalling swga.")
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


