#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h> // for mkdir
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>
#ifndef OSX
#include <sys/sysinfo.h> // to determine system memory
#endif
#include <sys/statvfs.h> // to determine available disk space
#include <dirent.h> // to clear the temp directory
#include <libgen.h> // for basename()
#include <string>

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

int max_memory; // the most memory dsk should alloc at any time, in MB
int max_disk_space; // the most disk space dsk should use at any time, in MB
extern bool output_histo;
#include "minia/Bank.h"
#include "minia/Pool.h"
#include "minia/Bloom.h"
#include "minia/Utils.h"
#include "minia/SortingCount.h"
#include "minia/Kmer.h"

int main(int argc, char *argv[])
{
    if(argc <  3)
    {
        fprintf (stderr,"%s: [d]isk [s]treaming of [k]-mers (constant-memory k-mer counting)\n",argv[0]);
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size [-t min_abundance] [-m max_memory] [-d max_disk_space] [-o out_prefix] [-histo] [-c] [-b]\n",argv[0]);
        fprintf (stderr,"Input file can be fasta, fastq, gzipped or not, or a text file containing one file name per line.\ndetails:\n [-t min_abundance] filters out k-mers seen ( < min_abundance ) times, default: 1 (all kmers are returned)\n [-m max_memory] is in MB, default: min(total system memory / 2, 5 GB) \n [-d max_disk_space] is in MB, default: min(available disk space / 2, reads file size)\n [-o out_prefix] saves results in [out_prefix].solid_kmers. default out_prefix = basename(input_file)\n [-histo] outputs histogram of kmers abundance\n [-c] write a Minia-compatible output file, i.e. discard k-mer counts\n [-b] use existing binary-converted reads, default: no (always recompute binary reads)\n");
#ifdef SVN_REV
fprintf(stderr,"Running dsk version %s\n",STR(SVN_REV));
#endif
        return 0;
    }

    // reads file
    Bank *Reads = new Bank(argv[1]);

    if (argv[2][0] == '-')
    {
        printf("please specify a k value\n");
        exit(1);
    }

    // kmer size
    sizeKmer = atoi(argv[2]);
    if (sizeKmer>(int)(sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %lu\n",sizeof(kmer_type)*4);
        exit(1);
    }

    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    // default solidity 
    nks = 1;

    // default max memory
    max_memory = 4*1024;
    #ifndef OSX
    struct sysinfo info;
    sysinfo(&info);
    int total_ram = (int)(((double)info.totalram*(double)info.mem_unit)/1024/1024);
    printf("Total RAM: %d MB\n",total_ram);
#else
    int total_ram = 128*1024;
#endif


    // default prefix is the reads file basename
    char *reads_path=strdup(argv[1]);
    string reads_name(basename(reads_path)); // posix basename() may alter reads_path
    free(reads_path);
    int lastindex = reads_name.find_last_of("."); 
    strcpy(prefix,reads_name.substr(0, lastindex).c_str()); 

    int verbose = 0;
    bool write_count = true;
    bool skip_binary_conversion = false;

    max_disk_space = 0;

    output_histo =false;
    for (int n_a = 3; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"-t")==0)
            nks = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-o")==0)
            strcpy(prefix,argv[n_a+1]);
 
        if (strcmp(argv[n_a],"-m")==0)
            max_memory = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-d")==0)
            max_disk_space = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-v")==0)
            verbose = 1;

        if (strcmp(argv[n_a],"-vv")==0)
            verbose = 2;
        
        if (strcmp(argv[n_a],"-histo")==0)
            output_histo = true;

        if (strcmp(argv[n_a],"-c")==0)
            write_count = false;
        
        if (strcmp(argv[n_a],"-b")==0)
            skip_binary_conversion = true;

        if (strcmp(argv[n_a],"--optimism")==0) // it's a shadow argument for now 
            optimism = atoi(argv[n_a+1]);
    }

    if (max_memory > total_ram)
    {
        printf("Maximum memory (%d MB), exceeds total RAM (%d MB). Setting maximum memory to %d MB.\n",max_memory,total_ram,total_ram/2);
        max_memory = total_ram/2;
    }

    STARTWALL(0);

    sorting_count(Reads, prefix, max_memory, max_disk_space, write_count, verbose, skip_binary_conversion);

    STOPWALL(0,"Total");

    delete Reads;

    return 0;
}


