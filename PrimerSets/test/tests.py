import tempfile
import unittest
from collections import namedtuple, OrderedDict
from StringIO import StringIO
import PrimerSets as ps
from PrimerSets.filter_primers import filter_primers
from PrimerSets.find_sets import find_sets

class ReadPrimersTests(unittest.TestCase):
    '''
    Testing primer read/writing functions
    '''
    def setUp(self):
        pass

    def testParsePrimer(self):
        primer = ps.parse_primer("""ATGC 100 200 0.5""")
        self.assertEqual(primer, ps.Primer(1, "ATGC", 200, 100, 0.5))

    def testReadPrimerFile(self):
        infile = StringIO("""AAAA 1 2 3
TTTT 1 2 3 4
GGGG 1 2 3""")
        primers = ps.read_primer_file(infile, False, True)
        self.assertEqual([ps.Primer(1, "AAAA", 2, 1, 3.0),
        ps.Primer(3, "GGGG", 2, 1, 3.0)], primers)


class FilterPrimersTests(unittest.TestCase):
    '''
    Testing filter_primer behavior
    '''
    def setUp(self):
        mockArgs = namedtuple("mockargs",
        ["max_bg_binding", "num_primers", "input", "output", "quiet"])
        self.tmpin = StringIO("""Some header file here
ATGC 50 99 0.5
TCGA 100 200 0.4
CGTA\t100\t98\t0.2
GCCT;1;2;3
CTTA 4 99 0.4
""")
        self.tmpout = tempfile.TemporaryFile()
        self.args = mockArgs(100, 2, self.tmpin, self.tmpout, True)


    def test_filter_primers(self):
        correct_result = [ps.Primer(4, "CGTA", 98, 100, 0.2),
        ps.Primer(6, "CTTA", 99, 4, 0.4)]
        test_result = filter_primers(self.args)
        self.assertEqual(correct_result, test_result)


    def tearDown(self):
        self.tmpout.close()


class GraphTests(unittest.TestCase):
    '''
    Testing mk_primer_graph and find_sets behavior
    '''
    def setUp(self):
        self.good_set = [ps.Primer(1, "ATGC", 1400, 15, 0),
                         ps.Primer(2, "GGCC", 1500, 10, 0),
                         ps.Primer(3, "CCTA", 2, 4, 0)]
        self.good_edges = [[1, 2], [1, 3]]
        self.good_result = '''p sp 3 2
n 1 1400
n 2 1500
n 3 2
e 1 2
e 1 3
'''
        self.bad_set = [(1,2,3,4), (4,3,2,1)]
        self.bad_edges = [[1], [1,2]]
        self.bad_handle = "tempfile!"
        self.tmp_infile = tempfile.NamedTemporaryFile()
        self.tmp_outfile = tempfile.NamedTemporaryFile()

    def test_correct_write(self):
        with tempfile.TemporaryFile() as tmp:
            ps.write_graph(self.good_set, self.good_edges, tmp)
            tmp.flush()
            tmp.seek(0)
            result = tmp.read()
            assert result == self.good_result

    def test_bad_handle(self):
        self.failUnlessRaises(ValueError, ps.write_graph, self.good_set,
        self.good_edges, self.bad_handle)

    def test_bad_primers(self):
        with tempfile.TemporaryFile() as tmp:
            self.failUnlessRaises(ValueError, ps.write_graph, self.bad_set,
            self.good_edges, tmp)

    def test_bad_edges(self):
        with tempfile.TemporaryFile() as tmp:
            self.failUnlessRaises(ValueError, ps.write_graph, self.good_set,
            self.bad_edges, tmp)

    def test_set_finder(self):
        '''
        Gives the set_finder a simple graph with two possible cliques, between
        edges {1,2,3} and {4,5,6}, where each edge weight is equal to its index
        and the background genome length is set to 10, and min_bg_bind_dist = 2.
        The only clique of the two that should satisfy those reqs is the first
        one, with a bg_ratio of 10/sum(1,2,3) = 1 (bc of integer division).
        '''
        test_graph = '''p sp 6 7
n 1 1
n 2 2
n 3 3
n 4 4
n 5 5
n 6 6
e 1 2
e 1 3
e 2 3
e 2 4
e 4 5
e 4 6
e 5 6
'''
        self.tmp_infile.write(test_graph)
        self.tmp_infile.flush()
        self.tmp_outfile = tempfile.NamedTemporaryFile()
        argvals = OrderedDict({'set_finder': 'set_finder',
            'min_bg_bind_dist': 2,
            'bg_genome_len': 10,
            'min_size': 3,
            'max_size': 3,
            'input': self.tmp_infile.name,
            'output': self.tmp_outfile.name})
        args = namedtuple('args', argvals.keys())
        argvals = args(**argvals)
        find_sets(argvals)
        result = self.tmp_outfile.read()
        self.assertEqual(result, "1,2,3 1\n")

    def tearDown(self):
        self.tmp_outfile.close()
        self.tmp_infile.close()


class HeterodimerTests(unittest.TestCase):
    '''
    Testing heterodimer checks
    '''
    def setUp(self):
        self.primers = []
        self.primers.append(ps.Primer(1, "ATGCTC", 0, 0, 0))
        # 4 contiguous
        self.primers.append(ps.Primer(2, "CAGCAT", 0, 0, 0))
        # 3 contiguous, alt
        self.primers.append(ps.Primer(3, "GAGGTA", 0, 0, 0))
        # 3 contiguous, alt 2
        self.primers.append(ps.Primer(4, "ATCGAC", 0, 0, 0))
        # valid
        self.primers.append(ps.Primer(5, "TTCCAC", 0, 0, 0))


    def test_compatible_primers(self):
        edges = ps.test_pairs(self.primers[0:2], 3)
        self.assertEqual(edges, [])
        edges = ps.test_pairs([self.primers[0], self.primers[4]], 3)
        self.assertEqual(edges, [[1, 5]])


class PrimerBindingLocationTests(unittest.TestCase):
    '''
    Testing functions to find primer binding locations
    '''
    def setUp(self):
        pass

    def test_find_locations(self):
        test_str = "TATATATAT"
        substr = "TAT"
        locs = ps.find_locations(substr, test_str)
        self.assertEqual(locs, [0, 2, 4, 6])

    def test_find_primer_locations(self):
        test_str = ">TATATATAT>TATA"
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write(test_str)
            tmp.flush()
            filename = tmp.name
            primer = ps.Primer(1, "TATA", 0, 0, 0)
            _, test_locs = ps.find_primer_locations(primer, filename)
        locs = sorted([0, 10, 1, 3, 5, 11, 2, 4, 6])
        self.assertEqual(test_locs, locs)


# class FgBindDistanceTest(unittest.TestCase):
#
#     def setUp(self):
#         self.line = "5 5 1 2 3"
#         self.locations = {1:{'seq':'',
#                              'loc':[1, 2, 3]},
#                           2:{'seq':'',
#                              'loc':[4, 5, 6]},
#                           3:{'seq':'',
#                              'loc':[7, 8, 10]}}
#
#
#     def test_fg_bind_dist(self):
#         primer_set, primers, max_dist, stdev = ps.fg_bind_distances(self.line, self.locations, ps.stdev)
#         self.assertEqual(max_dist, 2)
#         self.assertEqual(round(stdev, 4), 2.9345)
#         self.assertEqual(primer_set, [1, 2, 3])





# def test_test_pairs():
#     pass

# def test_max_consecutive_binding():
#     pass

if __name__ == "__main__":
    unittest.main()
