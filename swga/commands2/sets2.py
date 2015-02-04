import os
import subprocess
import swga
import swga.primers
import swga.graph as graph
from swga.commands2 import Command
from pkg_resources import resource_filename
from swga.primers import Primer, db

graph_fname = "compatibility_graph.dimacs"


def main(argv, cfg_file):
    cmd = Command('sets', cfg_file=cfg_file)
    cmd.parse_args(argv)
    print cmd.args, cfg_file
    find_sets(**cmd.args)


def find_sets(kmer_dir,
              output,
              min_size,
              max_size,
              max_hetdimer_bind,
              min_bg_bind_dist,
              bg_genome_len):
    assert kmer_dir
    primerdb_fp = os.path.join(kmer_dir, swga.primers.db_fname)
    if not os.path.isfile(primerdb_fp):
        swga.swga_error("Missing primers database: re-run `swga count`")

    db.init(primerdb_fp)
    
    # Reset all the primer IDs (only used for set_finder)
    Primer.update(pid = -1).execute()

    primers = list(Primer
                   .select()
                   .where(Primer.active==True)
                   .order_by(Primer.ratio.desc())
                   .execute())
    
    for i, primer in enumerate(primers):
        primer.pid = i+1
        primer.save()

    edges = graph.test_pairs(primers, max_hetdimer_bind)
    with open(graph_fname, 'wb') as out:
        graph.write_graph(primers, edges, out)
    
    set_finder = resource_filename("swga", "bin/set_finder")


    find_set_cmd = [set_finder, '-q', '-q', '-B', min_bg_bind_dist,
                    '-L', bg_genome_len, '-m', min_size, '-M',
                    max_size, '-a', '-u', '-r', 'unweighted-coloring',
                    graph_fname]
    find_set_cmd = " ".join([str(_) for _ in find_set_cmd])
    swga.message("> '%s'" % find_set_cmd)
    sf_proc = subprocess.Popen(args=find_set_cmd, shell=True, stdout=subprocess.PIPE,
                               stderr=None)
    for line in sf_proc.stdout.readlines():
        primerset, weight = swga.score.read_set_finder_line(line)
