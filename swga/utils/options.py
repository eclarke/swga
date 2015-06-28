import sys
import argparse
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
            prog='swga ' + name,
            description=description)
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

        # Setting some defaults
        opttype = argvals.get('opttype', 'opt')
        action = argvals.get('action', 'store')
        prefix = "--"
        help = argvals.get('desc', '')
        choices = argvals.get('choices', None)
        nargs = argvals.get('nargs', None)

        # Converting from 'std[in/out] string to actual streams
        default = argvals.get('default')
        if default == 'stdin':
            default = sys.stdin
        elif default == 'stdout':
            default = sys.stdout
        else:
            default = False

        # What kind of values can the argument take?
        _type = argvals.get('type')
        if _type == "File-w":
            _type = argparse.FileType('w')
        elif _type == "File-r":
            _type = argparse.FileType('r')
        elif _type:
            _type = eval(_type)

        # We remove the -- prefix for arguments (as opposed to options)
        if opttype == 'arg':
            prefix = ""

        # Flag-type options cannot have a type (even if None) specified in the
        # constructor
        if opttype == 'flag':
            parser.add_argument(
                prefix + argname,
                action='store_true',
                help=help)
        # We only use a non-False default value for stdin/stdout options- all
        # the other options use defaults from the config file for clarity
        elif default:
            parser.add_argument(
                prefix + argname,
                action=action,
                nargs=nargs,
                choices=choices,
                default=default,
                type=_type,
                help=help)
        # No default specified (annoyingly, specifying None conflicts with
        # specifying a type, so we must omit it from the constructor)
        else:
            parser.add_argument(
                prefix + argname,
                action=action,
                nargs=nargs,
                choices=choices,
                type=_type,
                help=help)

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
