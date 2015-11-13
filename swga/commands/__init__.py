from command import Command
from specfiles import command_names

map(__import__, command_names)
