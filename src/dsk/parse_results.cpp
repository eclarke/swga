#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <iostream>

#include "minia/Bank.h"
#include "minia/Kmer.h"
#include "minia/Utils.h"
#include "minia/SortingCount.h" // for uint_abundance_t

using namespace std;

char *solid_kmers_with_count_filename;
int coverage_threshold = 1;
bool dot_format = false;

struct kmer_count_raw_t
{
    unsigned char raw_data[sizeof(uint_abundance_t)+sizeof(kmer_type)];
};

int main(int argc, char *argv[])
{
    
    if(argc <  2)
    {
        fprintf (stderr,"parameters:\n");
        fprintf (stderr," dsk_solid_kmers_file [count threshold] [--dot]\n");
        return 0;
    }

    // first arg: kmer counts from DSK
    solid_kmers_with_count_filename = argv[1];

    // second arg: threshold 
    if (argc > 2)
    {
        if (argv[2][0]!='-')
        {
            coverage_threshold = atoi(argv[2]);
            if (coverage_threshold < 1)
            {
                printf("Coverage threshold must be > 0.\n");
                exit(1);
            }
        }
    }

    // optional argument: dot format
    for (int n_a = 2; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"--dot")==0)
            dot_format = true;
    }

    // get k from the solid k-mers file
    BinaryBank *SolidKmersWithCount = new BinaryBank(solid_kmers_with_count_filename, sizeof(uint_abundance_t) + sizeof(kmer_type),false);
    uint32_t kmer_nbits;
    SolidKmersWithCount->read(&kmer_nbits, 4);
    uint32_t k;
    SolidKmersWithCount->read(&k, 4);

    // set k-mer size now that we know it (and mask: Minia stuff)
    sizeKmer = k;
    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    if (kmer_nbits/2 != (sizeof(kmer_type)*4))
    {
       printf("The solid k-mers file was created with `make k=%d` but this software was compiled with `make k=%d`\n", kmer_nbits/2, sizeof(kmer_type)*4);
       exit(1);
    }

    char buffer[sizeKmer+1];
    kmer_count_raw_t kmer_count_raw;

    while (SolidKmersWithCount->read_element_buffered(&(kmer_count_raw.raw_data[0])))
    {
        kmer_type *kmer = (kmer_type*)&(kmer_count_raw.raw_data[0]);
        uint_abundance_t *abundance = (uint_abundance_t*) &(kmer_count_raw.raw_data[sizeof(kmer_type)]);
        if (*abundance < coverage_threshold)
            continue;
        code2seq(*kmer, buffer, sizeKmer, kmerMask);

        if (dot_format)
        {
            string lower = buffer;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            cout << lower << ";\n";
        }
        else
            cout << buffer << " " << *abundance << "\n";
    }
    SolidKmersWithCount->close();

    return 0;
}


