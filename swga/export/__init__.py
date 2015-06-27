import os

from bedfile import BedFile
from bedgraph import BedGraph


def _mk_folder(folder, fg_genome_fp, setstr):
    fg_name = ".".join(
        os.path.basename(fg_genome_fp).split(".")[:-1]) + "_export"
    output_folder = os.path.join(folder, fg_name, setstr)
    if not os.path.isdir(output_folder):
            os.makedirs(output_folder)
    return output_folder
