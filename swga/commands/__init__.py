import argutils
import argutils.read
import argutils.export
import swga.utils as utils
import activate
import count
import export
import filter
import find_sets
import score
import summary


_commands = [
    'activate',
    'count',
    'export',
    'filter',
    'find_sets',
    'score',
    'summary',
]


def create_config_file():
    cfg_file_str = ""
    for cmd in _commands:
        spec = utils.specfile(cmd)
        opts = argutils.read.from_yaml(spec)
        if opts:
            if opts['_meta'].get("_exclude"):
                continue
            cfg_file_str += argutils.export.to_config(cmd, opts) + "\n"
    return cfg_file_str


Activate = activate.Activate
Count = count.Count
Export = export.Export
Filter = filter.Filter
FindSets = find_sets.FindSets
Score = score.Score
Summary = summary.Summary
