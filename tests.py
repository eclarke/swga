from primer_cliques import *

# test read_primers
primers = read_primers('test_primers_binding.txt')
assert len(primers) == 200
assert primers[6].id == 6
assert primers[6].seq == "GTGTGTG"
assert primers[6].bg_freq == 1093426
