import subprocess
from swga.commands import Command
from pkg_resources import resource_filename


def main(argv, cfg_file):
    cmd = Command('sets', cfg_file = cfg_file)
    cmd.parse_args(argv)
    find_sets(**cmd.args)


def find_sets(input,
              output,
              min_size,
              max_size,
              min_bg_bind_dist,
              bg_genome_len):
    '''
    Calls the set_finder binary with the specified options on the
    heterodimer compatibility graph and outputs valid sets for
    post-processing.
    '''
    set_finder = resource_filename("swga", "bin/set_finder")

    output = '> ' + output if output else ''
    find_set_cmd = [set_finder, '-q', '-q', '-B', min_bg_bind_dist,
                    '-L', bg_genome_len, '-m', min_size, '-M',
                    max_size, '-a', '-u', '-r', 'unweighted-coloring',
                    input, output]
    find_set_cmd = [str(_) for _ in find_set_cmd]

    subprocess.check_call(" ".join(find_set_cmd), shell=True)
