import os
import yaml
from pkg_resources import resource_stream
from collections import OrderedDict
import argutils

command_names = [
	'activate',
	'count',
	'export',
	'filter',
	'find_sets',
	'init',
	'score',
	'summary'
]

def get_cmd_specfile(cmdname):
	fp = os.path.join('commands', 'specfiles', cmdname+'.yaml')
	return resource_stream("swga", fp)