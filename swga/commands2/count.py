import os
import subprocess
import swga
from pkg_resources import resource_filename
from swga.commands2 import Command
import swga.resources as resources

def main(argv, cfg_file):
    cmd = Command('count', cfg_file = cfg_file)
    cmd.parse_args(argv)
    count_mers(**cmd.args)


def count_mers(fg_genome_fp, bg_genome_fp, min_size, max_size, threshold, output_dir):
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
    cmdstr = "{dsk} {{genome}} {{i}} -o {{output}}".format(dsk=dsk)
    parser_cmdstr = ("{parse_dsk} {{output}}.solid_kmers_binary "
                     "{threshold} > {{output}}").format(parse_dsk=parse_dsk, 
                                                        threshold=threshold)
    # Create the directories if they don't exist already
    for d in (output_dir, fg_kmer_dir, bg_kmer_dir):
        swga.mkdirp(d)

    tmp_files = []
    for i in xrange(min_size, max_size+1):
        output_name = '%i-mers.fa' % i
        fg_cmdstr = cmdstr.format(i = i, genome = fg_genome_fp, 
                                  output = output_name)
        fg_parser_cmdstr = parser_cmdstr.format(output = output_name)

        bg_cmdstr = cmdstr.format(i = i, genome = bg_genome_fp,
                                  output = output_name)
        bg_parser_cmdstr = parser_cmdstr.format(output = output_name)

        print fg_cmdstr
        subprocess.check_call(fg_cmdstr, shell=True, cwd=fg_kmer_dir)
        tmp_files.append(os.path.join(fg_kmer_dir, output_name+'.solid_kmers_binary'))
        print fg_parser_cmdstr
        subprocess.check_call(fg_parser_cmdstr, shell=True, cwd=fg_kmer_dir)
        print bg_cmdstr
        subprocess.check_call(bg_cmdstr, shell=True, cwd=bg_kmer_dir)
        tmp_files.append(os.path.join(bg_kmer_dir, output_name+'.solid_kmers_binary'))
        print bg_parser_cmdstr
        subprocess.check_call(bg_parser_cmdstr, shell=True, cwd=bg_kmer_dir)
        
    subprocess.check_call(["rm *.solid_kmers_binary"], shell=True, cwd=fg_kmer_dir)
    subprocess.check_call(["cat *-mers.fa > fg_mers.fa"], shell=True, cwd=fg_kmer_dir)
    subprocess.check_call(["rm *.solid_kmers_binary"], shell=True, cwd=bg_kmer_dir)
    subprocess.check_call(["cat *-mers.fa bg_mers.fa"], shell=True, cwd=bg_kmer_dir)

if __name__ == '__main__':
    main(None, None)
