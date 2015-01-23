from pkg_resources import resource_stream
import yaml
from collections import OrderedDict

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_constructor(_mapping_tag, dict_constructor)
yaml.add_representer(OrderedDict, dict_representer)

def load_swga_opts():
    with resource_stream("swga", "data/options.yaml") as opts_fp:
        opts = yaml.load(opts_fp)
        return opts


def opts2cfg(opts):
    out_str     = ""
    section_str = "{desc}\n[{section}]\n"
    option_str  = "# {desc}\n{opt} = {default}\n"
    default_str = "{item} = {default}\n"

    for section in opts.keys():
        desc = opts[section].get("desc")
        desc = "\n## "+desc if desc else ""
        out_str += section_str.format(desc=desc, section=section)

        if section == "DEFAULT":
            for item in opts[section]:
                out_str += default_str.format(item=item,
                                              default=opts[section][item])
            continue
        
        for opt in opts[section].keys():
            if opt == "desc": continue
            option = opts[section][opt]
            desc = option.get("desc")
            default = option.get("default")
            out_str += option_str.format(desc=desc, opt=opt, default=default)
    
    return(out_str)
            
