import subprocess
import swga
import swga.primers
import swga.graph as graph
from swga.commands import Command
from pkg_resources import resource_filename
from swga.database import Primer, update_in_chunks

graph_fname = "compatibility_graph.dimacs"


def main(argv, cfg_file):
    cmd = Command('find_sets', cfg_file=cfg_file)
    cmd.parse_args(argv)
    swga.primers.init_db(cmd.primer_db)
    find_sets(cmd.min_bg_bind_dist, cmd.bg_genome_fp,
              cmd.min_size, cmd.max_size, cmd.bg_genome_len)
    make_graph(cmd.max_hetdimer_bind, graph_fname)

    
def make_graph(max_hetdimer_bind, outfile):
    '''Selects all active primers and outputs a primer compatibility graph.'''

    # Reset all the primer IDs (ids are only used for set_finder)
    Primer.update(pid = -1).execute()

    primers = (Primer.select().select(Primer.active==True)
               .order_by(Primer.ratio.desc()))
    
    if len(primers) == 0:
        swga.swga_error("No active sets found. Run `swga filter` first.")

    for i, p in enumerate(primers):
        p.pid = i + 1

    update_in_chunks(primers, show_progress=False)

    swga.message("Composing primer compatibility graph...")
    edges = graph.test_pairs(primers, max_hetdimer_bind)

    if len(edges) == 0:
        swga.swga_error("No compatible primers. Try relaxing your parameters.")

    with open(outfile, 'wb') as out:
        graph.write_graph(primers, edges, out)

        
def find_sets(min_bg_bind_dist, bg_genome_fp, min_size, max_size, bg_genome_len):
    swga.message("Now finding sets. If nothing appears, try relaxing your parameters.")
    set_finder = resource_filename("swga", "bin/set_finder")
    find_set_cmd = [set_finder, '-q', '-q', '-B', min_bg_bind_dist,
                    '-L', bg_genome_len, '-m', min_size, '-M',
                    max_size, '-a', '-u', '-r', 'unweighted-coloring',
                    graph_fname]
    find_set_cmd = " ".join([str(_) for _ in find_set_cmd])
    subprocess.check_call(find_set_cmd, shell=True)
