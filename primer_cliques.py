## primer_cliques.py 
# Uses cliquer to find maximal cliques (completely connected subgraphs) of 
# primers. 
#
# Erik Clarke - ecl@mail.med.upenn.edu
from itertools import combinations
from argparse import ArgumentParser
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
import sys


def write_graph(num_nodes, arcs, fname):
    '''
    Writes graph in DIMACS graph format, specified in the cliquer user manual.
    See http://users.tkk.fi/~pat/cliquer.html
    
    An arc is simply a list of the form [first_node, second_node]
    'arcs' is a list of arcs
    '''
    with open(fname, 'w') as output:
        output.write('p sp {} {}\n'.format(num_nodes, len(arcs)))
        for arc in arcs:
            output.write('e {} {} \n'.format(arc[0], arc[1]))


def read_primers(primer_fp, num_primers, bg=None):
    '''
    Reads first <num_primers> lines of a tab-delimited file where the first 
    column contains the primer sequences.
    If bg != None, checks background binding frequencies.
    Returns a list of primers in format [[primer_id primer_seq ]...]
    '''
    primers = []
    pool = None
    with open(primer_fp) as primerfile:
        for i, line in enumerate(primerfile.readlines()):
            if i >= num_primers: break
            primer = line.strip('\n').split('\t')[0]
            primers.append((i, primer))
    if bg:
        print ('\t Determining background binding counts...')
        bg_frequencies = [bg_binding_freq(primer, bg) for primer in primers]
    return(primers, bg_frequencies)


def bg_binding_freq_fp(primer, bg_fp):
    '''
    Searches for the primer in the background genome, specified by the
    filename in bg_fp, using a sliding window of the size of the primer.
    '''
    count = 0
    psize = len(primer) # assuming len of primer == number of bytes
    i = 0
    with open(bg_fp) as bg:
        bg.readline()
        start = bg.tell() # what's the index after the first line?
        while True:
            chunk = bg.read(psize)
            # if the size is smaller than we requested, we've hit the end of
            # the file
            if len(chunk) < psize:
                return count
            if chunk[-1] == '\n':
                chunk = chunk + bg.read(1)
                chunk = chunk.replace('\n', '')
            if chunk == primer:
                count += 1
            i += 1
            bg.seek(start+i)


def bg_binding_freq(primer, bg):
    count = start = 0
    while True:
        start = bg.find(primer, start) + 1
        if start > 0:
            count += 1
        else:
            return count


def test_pairs(starting_primers, max_binding):
    '''
    Adds a primer pair to the list of "arcs" if it passes the heterodimer
    filter using the max_binding cutoff.
    '''
    arcs = []
    for x, y in combinations(starting_primers, 2):
        if max_consecutive_binding(x[1], y[1]) < max_binding:
            arcs.append([x[0], y[0]])
    return arcs


def max_consecutive_binding(mer1, mer2):
    '''
    Return the maximum number of consecutively binding mers
    when comparing two different mers, using the reverse compliment.
    '''
    binding = { 'A': 'T', 'T': 'A',
                'C': 'G', 'G': 'C',   
                '_':  False}

    # Swap variables if the second is longer than the first
    if len(mer2) > len(mer1):
        mer1, mer2 = mer2, mer1

    # save the len because it'll change when we do a ljust
    mer1_len = len(mer1)
    # reverse mer2,
    mer2 = mer2[::-1]
    # pad mer one to avoid errors
    mer1 = mer1.ljust(mer1_len + len(mer2), "_")

    max_bind = 0
    for offset in range(mer1_len):
        consecutive = 0
        for x in range(len(mer2)):
            if binding[mer1[offset+x]] == mer2[x]:
                consecutive += 1
                if consecutive > max_bind:
                    max_bind = consecutive
            else:
                consecutive = 0
    return max_bind


def main():
    #--- Argument parsing ---#
    description_string = '''Reads in primers from a tab-delimited file, performs heterodimer checks and calculates background genome binding frequency, then writes results in DIMACS format for use with cliquer.c.'''
    parser = ArgumentParser(description = description_string)
    parser.add_argument('-n', '--num_primers', help="Desired number of primers (read from top of input file", type=int, default=50)
    parser.add_argument('-m', '--max_binding', help="Max number of consecutive complimentary bases allowed in the heterodimer filter", type=int, default=3)
    parser.add_argument('-o', '--output', help="Filename to store the DIMACS output graph.", required=True)
    parser.add_argument('-b', '--bg_genome', help="Location of the FASTA-format background genome. We'll be reading this all into memory. If unspecified, skip background checking.", default=None)
    parser.add_argument('primers', action="store")
    args = vars(parser.parse_args())


    
    num_primers = args["num_primers"]
    max_binding = args["max_binding"]
    bg_fp = args["bg_genome"]
    bg = None
    if bg_fp:
        print("Reading in background genome...")
        with open(bg_fp) as bg_handle:
            bg = "".join(_.strip('\n') for i, _ in enumerate(bg_handle.readlines()) if i > 0)
        print("\tdone.")
    print("Reading in primers and assigning bg hit counts...")
    primers, bg_binding_counts = read_primers(args['primers'],
                                              num_primers, bg)
    print("\tdone.")
    print("Testing all primer pairs in heterodimer filter...")
    arcs = test_pairs(primers, max_binding)
    print("\tdone.")
    write_graph(num_primers+1, arcs, 'test_graph.gr')
    print("Wrote graph to test_graph.gr.")
    print("Done!")

if __name__ == '__main__':
        main()




