import sys,struct

if len(sys.argv) < 2:
    sys.exit( "usage: [solid_kmers_binary file]" )

fin = open(sys.argv[1], "rb")
kmer_nbits = struct.unpack("i",fin.read(4))[0]
k = struct.unpack("i",fin.read(4))[0]

try:
    while True:
        kmer_binary = struct.unpack('B'*(kmer_nbits / 8), fin.read( kmer_nbits / 8) )
        abundance = struct.unpack( 'I', fin.read(4) ) [0]
        kmer = ""
        for i in xrange(k):
            kmer = "ACTG"[(kmer_binary[i/4] >> (2*(i%4)) ) % 4] + kmer
        print kmer, abundance
except:
    fin.close()
