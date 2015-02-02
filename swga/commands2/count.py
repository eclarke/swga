#from __future__ import division
import os
import subprocess
import swga
import csv
import struct
from swga.commands2 import Command
import swga.resources as resources
import sqlite3
from swga.primers import write_primer_file, mk_primer_tbl, update_primer_tbl, Primer


def main(argv, cfg_file):
    cmd = Command('count', cfg_file=cfg_file)
    cmd.parse_args(argv)
    count_mers(**cmd.args)


def count_mers(fg_genome_fp,
               bg_genome_fp, 
               min_size, 
               max_size, 
               threshold, 
               output_dir):
    assert os.path.isfile(fg_genome_fp)
    assert os.path.isfile(bg_genome_fp)
    try:
        dsk, parse_dsk = resources.get_dsk()
    except ValueError as e:
        swga.swga_error("Error: %s \n\n Required dsk binaries not found. "
                        "Try re-installing SWGA." % e.message)

    fg_genome_name = fg_genome_fp.split(os.sep).pop()
    bg_genome_name = bg_genome_fp.split(os.sep).pop()
    fg_kmer_dir = os.path.join(output_dir, fg_genome_name+'.fgkmers')
    bg_kmer_dir = os.path.join(output_dir, bg_genome_name+'.bgkmers')
    cmdstr = "{dsk} {{genome}} {{i}} -o {{output}} -t {threshold}".format(
        dsk=dsk, threshold=threshold)
    parser_cmdstr = ("{parse_dsk} {{output}}.solid_kmers_binary "
                     "> {{output}}").format(parse_dsk=parse_dsk, 
                                            threshold=threshold)
    # Create the directories if they don't exist already
    for d in (output_dir, fg_kmer_dir, bg_kmer_dir):
        swga.mkdirp(d)

    tmp_files = []
    for i in xrange(min_size, max_size+1):
        output_name = '%i-mers.fa' % i

        fg_cmdstr = cmdstr.format(i=i, genome=fg_genome_fp,
                                  output=output_name)
        fg_parser_cmdstr = parser_cmdstr.format(output=output_name)

        bg_cmdstr = cmdstr.format(i=i, genome=bg_genome_fp,
                                  output=output_name)
        bg_parser_cmdstr = parser_cmdstr.format(output=output_name)

        # Runs each command and adds the resulting file to the tmp_files array
        tmp_files.append(
            _run_and_log(fg_cmdstr, output_name+'.solid_kmers_binary',
                         fg_kmer_dir))
        tmp_files.append(
            _run_and_log(fg_parser_cmdstr, output_name, fg_kmer_dir))

        tmp_files.append(
            _run_and_log(bg_cmdstr, output_name+'.solid_kmers_binary',
                         bg_kmer_dir))
        tmp_files.append(
            _run_and_log(bg_parser_cmdstr, output_name, bg_kmer_dir))

    subprocess.check_call(["cat *-mers.fa > fg_mers.fa"], shell=True, cwd=fg_kmer_dir)   
    subprocess.check_call(["cat *-mers.fa > bg_mers.fa"], shell=True, cwd=bg_kmer_dir)

    [os.remove(f) for f in tmp_files]

    swga.message("Merging primers...")
    primers = merge_mers(os.path.join(fg_kmer_dir, "fg_mers.fa"),
                         os.path.join(bg_kmer_dir, "bg_mers.fa"))

    with open(os.path.join(os.path.abspath(output_dir), 'primers.txt'), 'wb') as out:
        write_primer_file(primers, out, header=True)

    conn = sqlite3.connect("primers.sqlitedb")
    mk_primer_tbl(conn)
    update_primer_tbl(primers, conn)
    

def merge_mers(fg_mers_fp, bg_mers_fp):
    fg = bg = {}

    with open(fg_mers_fp) as _fg:
        rows = csv.DictReader(_fg, fieldnames=["seq", "freq"], delimiter=' ')
        fg = dict((row['seq'], row['freq']) for row in rows)
        
    with open(bg_mers_fp) as _bg:        
        rows = csv.DictReader(_bg, fieldnames=['seq', 'freq'], delimiter=' ')
        bg = dict((row['seq'], row['freq']) for row in rows)

    # Keep only the primers in fg
    primers = [Primer(id=i, seq=seq, fg_freq=fg[seq], bg_freq=bg.get(seq, 0))
               for i, seq in enumerate(fg.keys())]
    return primers


def parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        kmer_nbits = struct.unpack('i', f.read(4))[0]
        k = struct.unpack('i', f.read(4))[0]
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits / 8),
                                            f.read(kmer_nbits / 8))
                freq = struct.unpack('I', f.read(4))[0]
                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[(kmer_binary[i/4] >> (2 * (i%4))) % 4] + kmer
                yield kmer, freq
        except struct.error:
            pass
    
        
        
def _run_and_log(cmdstr, fname, cwd):
    swga.message(cmdstr)
    subprocess.check_call(cmdstr, shell=True, cwd=cwd)
    return os.path.abspath(os.path.join(cwd, fname))
