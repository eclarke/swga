import os
import subprocess
import swga
from pkg_resources import resource_filename
from swga.commands2 import Command


def main(argv, cfg_file):
    print "Count command"
    cmd = Command('count', cfg_file = cfg_file)
    cmd.parse_args(argv)
    print cmd.args
    count_mers(**cmd.args)


def count_mers(fg_genome_fp, bg_genome_fp, min_size, max_size, threshold, output_dir):
    dsk = resource_filename('swga', 'bin/dsk')
    parse_dsk = resource_filename('swga', 'bin/parse_results')
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

    for i in xrange(min_size, max_size+1):
        fg_output_name = bg_output_name = '%i-mers.fa' % i
        fg_cmdstr = cmdstr.format(i = i, genome = fg_genome_fp, 
                                  output = fg_output_name)
        fg_parser_cmdstr = parser_cmdstr.format(output = fg_output_name)

        bg_cmdstr = cmdstr.format(i = i, genome = bg_genome_fp,
                                  output = bg_output_name)
        bg_parser_cmdstr = parser_cmdstr.format(output = bg_output_name)

        print fg_cmdstr
        subprocess.check_call(fg_cmdstr, shell=True, cwd=fg_kmer_dir)
        print fg_parser_cmdstr
        subprocess.check_call(fg_parser_cmdstr, shell=True, cwd=fg_kmer_dir)
        print bg_cmdstr
        subprocess.check_call(bg_cmdstr, shell=True, cwd=bg_kmer_dir)
        print bg_parser_cmdstr
        subprocess.check_call(bg_parser_cmdstr, shell=True, cwd=bg_kmer_dir)
    
    subprocess.check_call(["rm *.solid_kmers_binary"], shell=True, cwd=fg_kmer_dir)
    subprocess.check_call(["cat *-mers.fa > fg_mers.fa"], shell=True, cwd=fg_kmer_dir)
    subprocess.check_call(["rm *.solid_kmers_binary"], shell=True, cwd=bg_kmer_dir)
    subprocess.check_call(["cat *-mers.fa bg_mers.fa"], shell=True, cwd=bg_kmer_dir)

if __name__ == '__main__':
    main(None, None)
