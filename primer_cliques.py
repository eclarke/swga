## primer_cliques.py 
# Uses cliquer to find maximal cliques (completely connected subgraphs) of 
# primers. 
#
# Erik Clarke - ecl@mail.med.upenn.edu
from itertools import combinations


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


def read_primers(primer_fp, num_primers, bg):
	'''
	Reads first <num_primers> lines of a tab-delimited file where the first 
	column contains the primer sequences.
	Returns a list of primers in format [[primer_id primer_seq]...]
	'''
	primers = []
	with open(primer_fp) as primerfile:
		for i, line in enumerate(primerfile.readlines()):
			if i > num_primers: break

			primers.append([i+1, line.strip('\n').split('\t')[0]])
	return(primers)


def test_pairs(starting_primers, max_binding):
	'''
	Adds a primer pair to the list of "arcs" if it passes the heterodimer
	filter using the max_binding cutoff.
	'''
	arcs = []
	for x, y in combinations(starting_primers, 2):
		if max_consecutive_binding(x[1], y[1]) < max_binding:
			arcs.append([x[0], y[0]])
		# else:
		# 	print("Primers {} and {} failed.".format(x[1], x[1]))
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


# def bg_binding_occurances(primer, bg):
# 	'''
# 	Waaaaayyyy faster than regex, but needs to load the bg genome in memory.
# 	'''
# 	count = start = 0
#     while True:
#         start = string.find(sub, start) + 1
#         if start > 0:
#             count+=1
#         else:
#             return count

def main():
	num_primers = 200
	max_binding = 3
	bg_fp = 'humangenome.fasta'
	bg = ""
	print("Reading in background genome...")
	# with open(bg_fp) as bg_handle:
	# 	bg = "".join(_strip('\n') for i, _ in enumerate(bg_handle.readlines()) if i > 0)
	print("\tdone.")
	print("Reading in primers and assigning bg hit counts...")
	primers = read_primers('selected-mers', num_primers, bg)
	print("\tdone.")
	print("Testing all primer pairs in heterodimer filter...")
	arcs = test_pairs(primers, max_binding)
	print("\tdone.")
	write_graph(num_primers+1, arcs, 'test_graph.gr')
	print("Wrote graph to test_graph.gr.")
	print("Done!")

if __name__ == '__main__':
	main()




