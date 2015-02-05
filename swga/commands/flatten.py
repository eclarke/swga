import subprocess
import sys
import os
from swga.commands import Command
    
        
def main(argv, cfg_file):
    cmd = Command('flatten', cfg_file = cfg_file)
    cmd.parse_args(argv)
    flatten(**cmd.args)


def flatten(input, output, force):
    if not os.path.isfile(input):
        sys.stderr.write("Error: input file specified does not "
                         "exist.\n")
        exit(1)
                         
    if not output:
        output = input+'.flattened'

    if os.path.isfile(output):
        if force:
            os.remove(output)
        else:
            sys.stderr.write("Error: output file exists.\n")
            exit(1)

    cmdstr = "sed 's/>.*/>/' {input} | tr -d '\\n' | tr '[:lower:]' '[:upper:]' > {output}"
    cmdstr = cmdstr.format(input=input, output=output)
    subprocess.check_call(cmdstr, shell=True)
    sys.stderr.write("{input} --> {output}\n".format(input=input, output=output))
    
