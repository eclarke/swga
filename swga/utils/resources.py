'''
Provides methods for retrieving SWGA resources
'''
import yaml
import swga
from collections import OrderedDict
from pkg_resources import (
    resource_exists, resource_filename, resource_stream
)


def _get_resource_file(fp):
    resource_path = resource_filename("swga", fp)
    if not resource_exists("swga", fp):
        swga.swga_error("Resource does not exist: {}".format(resource_path))    
    return resource_path


def get_dsk():
    return _get_resource_file("bin/dsk")


def get_setfinder():
    return _get_resource_file("bin/set_finder")


def get_swga_opts():
    with resource_stream("swga", "data/options.yaml") as opts_stream:
        opts = yaml.load(opts_stream)
        if opts == {}:
            swga.swga_error("Empty options.yaml file- reinstall SWGA.")
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


