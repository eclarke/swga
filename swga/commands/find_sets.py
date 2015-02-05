import os
import subprocess
import swga
import swga.primers
import swga.graph as graph
from swga.core import chunk_iterator
from swga.commands import Command
from pkg_resources import resource_filename
from swga.primers import Primer, upsert_chunk

graph_fname = "compatibility_graph.dimacs"


def main(argv, cfg_file):
    cmd = Command('find_sets', cfg_file=cfg_file)
    cmd.parse_args(argv)
    find_sets(**cmd.args)


def find_sets(primer_db,
              min_size,
              max_size,
              max_hetdimer_bind,
              min_bg_bind_dist,
              bg_genome_len):

    swga.primers.init_db(primer_db)
    
    # Reset all the primer IDs (only used for set_finder)
    Primer.update(pid = -1).execute()

    primers = list(Primer
                   .select()
                   .where(Primer.active==True)
                   .order_by(Primer.ratio.desc())
                   .execute())
    
    if len(primers) == 0:
        swga.swga_error("No active sets found. Run `swga filter` first.")

    for i, p in enumerate(primers):
        p.pid = i + 1
    chunk_iterator(primers, upsert_chunk, show_progress=False)
    # for i, primer in enumerate(primers):
    #     primer.pid = i+1
    #     primer.save()

    swga.message("Composing primer compatibility graph...")
    edges = graph.test_pairs(primers, max_hetdimer_bind)
    if len(edges) == 0:
        swga.swga_error("No compatible primers. Try relaxing your parameters.")
    with open(graph_fname, 'wb') as out:
        graph.write_graph(primers, edges, out)
    
    swga.message("Now finding sets. If nothing appears, try relaxing your parameters.")
    set_finder = resource_filename("swga", "bin/set_finder")
    find_set_cmd = [set_finder, '-q', '-q', '-B', min_bg_bind_dist,
                    '-L', bg_genome_len, '-m', min_size, '-M',
                    max_size, '-a', '-u', '-r', 'unweighted-coloring',
                    graph_fname]
    find_set_cmd = " ".join([str(_) for _ in find_set_cmd])
    swga.message("Set finder command:")
    swga.message(find_set_cmd)
    subprocess.check_call(find_set_cmd, shell=True)
