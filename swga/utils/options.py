import yaml
import sys
import argparse
from pkg_resources import resource_stream
from collections import OrderedDict
from swga.clint.textui import puts, max_width, indent
from StringIO import StringIO


def load_swga_opts():
    '''Returns the parsed options.yaml file as a Python object.'''
    with resource_stream("swga", "data/options.yaml") as opts_fp:
        opts = yaml.load(opts_fp)
        return opts


def cfg_from_opts(opts):
    '''
    Takes a dictionary of options (with desc and default fields) and creates a
    config file from them, parseable by ConfigParser.
    
    Returns: the contents of the config file, as a string.
    '''
    out_str     = ""
    section_str = "[{section}]\n"
    opt_str  = "{opt} = {default}\n"

    for section in opts.keys():
        desc = opts[section].get("desc")
        desc = "\n" + format_comment(desc, quote='##')
        out_str += desc + section_str.format(section=section)

        if section == "DEFAULT":
            for opt in opts[section]:
                out_str += opt_str.format(opt=opt, default=opts[section][opt])
            continue
        
        for opt in opts[section].keys():
            if opt == "desc": continue  # already handled this in above
            option = opts[section][opt]
            desc = format_comment(option.get("desc"))
            default = option.get("default")            
            out_str += desc + opt_str.format(opt=opt, default=default)
    
    return(out_str)


def mk_argparser_from_opts(opts, cmd_name, description=None):
    '''
    Uses a defined section (`cmd_name`) of the options to create an
    ArgumentParser object.

    Ignores the default value specified in opts (in preference to the defaults
    given by the config file later).

    Returns: the generated ArgumentParser
    '''

    def make_parser(name, description=None):
        description = description if description else ""
        parser = argparse.ArgumentParser(
            prog = 'swga ' + name,
            description = description)
        return parser

    if cmd_name not in opts:
        return make_parser(cmd_name, description)

    section = opts[cmd_name]

    # Use passed description if provided, else use the one from the opts
    description = description if description else section.get("desc")
    parser = make_parser(cmd_name, description)

    for argname in section:
        if argname == "desc": continue
        argvals = section[argname]
        _type = eval(argvals.get('type', 'str'))

        parser.add_argument(
            '--' + argname,
            action = argvals.get('action', 'store'),
            type = _type,
            help = argvals.get('desc', ''),
            required = argvals.get('required', False))

    return parser

def format_comment(desc, width=72, quote='#'):
    '''
    Wraps a description at the specified width and puts a comment character
    before each line. Returns a formatted string.
    '''
    if not desc:
        return ""
    string = StringIO()
    stream = string.write
    desc = max_width(desc, width)
    with indent(len(quote)+1, quote):
        puts(desc, stream=stream)
    return string.getvalue()


## This ensures that the options we get back from the yaml file are in order by
## replacing the default constructor with an OrderedDict
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def _dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def _dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_constructor(_mapping_tag, _dict_constructor)
yaml.add_representer(OrderedDict, _dict_representer)
