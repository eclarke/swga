import sys
import yaml
import argparse
from collections import OrderedDict
from swga.clint.textui import puts, max_width, indent
from StringIO import StringIO


def cfg_from_opts(opts):
    '''
    Takes a dictionary of options (with desc and default fields) and creates a
    config file from them, parseable by ConfigParser.
    
    Returns: the contents of the config file, as a string.
    '''
    out_str = ""
    section_str = "[{section}]\n"
    opt_str = "{opt} = {default}\n"

    for section in opts.keys():
        if section == "INTERNAL":
            continue
        meta = opts[section]["META"]
        if not meta.get('incfg', True):
            continue
        desc = meta.get('desc', '')
        desc = "\n" + _format_comment(desc, quote='##')
        out_str += desc + section_str.format(section=section)

        for opt in opts[section].keys():
            if opt == "META":
                continue
            if opt == "desc": 
                continue
            if opt == "incfg":
                continue    
            option = opts[section][opt]
            if not option.get("incfg", True):
                continue
            desc = _format_comment(option.get("desc"))
            default = option.get("default")
            if default == "None" or default is None:
                default = ""
            out_str += desc + opt_str.format(opt=opt, default=default)
    
    return(out_str)


def argparser_from_opts(opts, cmd_name, description=None):
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
    meta = opts[cmd_name]['META']
    # Use passed description if provided, else use the one from the opts
    description = description if description else meta['desc']
    parser = make_parser(cmd_name, description)

    for argname in section:

        if argname == "META":
            continue
        argvals = section[argname]
        _type = argvals.get('type')
        if _type == "File-w":
            _type = argparse.FileType('w')
        elif _type == "File-r":
            _type = argparse.FileType('r')
        elif _type:
            _type = eval(_type)

        opttype = argvals.get('opttype', 'opt')
        action = argvals.get('action', 'store')
        prefix = "--"
        help = argvals.get('desc', '')
        choices = argvals.get('choices', None)
        nargs = argvals.get('nargs', None)
        default = argvals.get('default')
        if default == 'stdin':
            default = sys.stdin
        elif default == 'stdout':
            default = sys.stdout
        else:
            default = False
        if opttype == 'arg':
            prefix = ""
        elif opttype == 'flag':
            action = 'store_true'
        if opttype == 'flag':
            print argname, default
            parser.add_argument(
                prefix + argname,
                action = action,
                help = help)
        elif default:
            print argname, default
            parser.add_argument(
                prefix + argname,
                action = action,
                nargs = nargs,
                choices = choices,
                default = default,
                type = _type,
                help = help)
        else:
            parser.add_argument(
                prefix + argname,
                action = action,
                nargs = nargs,
                choices = choices,
                type = _type,
                help = help)


    return parser


def _format_comment(desc, width=72, quote='#'):
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


# This ensures that the options we get back from the yaml file are in order by
# replacing the default constructor with an OrderedDict
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG


def _dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())


def _dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_constructor(_mapping_tag, _dict_constructor)
yaml.add_representer(OrderedDict, _dict_representer)
