import primer_sets as ps
from primer_sets import Primer
import tempfile
import unittest
from StringIO import StringIO

class GraphTests(unittest.TestCase):

    def setUp(self):
        self.good_set = [Primer(1, "ATGC", 1400, 15),
                         Primer(2, "GGCC", 1500, 10),
                         Primer(3, "CCTA", 2, 4)]
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
        
    def test_correct_write(self):
        with tempfile.TemporaryFile() as tmp:
            ps.write_graph(self.good_set, self.good_edges, tmp)
            tmp.flush()
            tmp.seek(0)
            result = tmp.read()
            assert result == self.good_result

    def test_bad_handle(self):
        self.failUnlessRaises(ValueError, ps.write_graph,
                              self.good_set, self.good_edges,
                              self.bad_handle) 

    def test_bad_primers(self):
        with tempfile.TemporaryFile() as tmp:
            self.failUnlessRaises(ValueError, ps.write_graph,
                                  self.bad_set, self.good_edges, tmp)                                   

    def test_bad_edges(self):
        with tempfile.TemporaryFile() as tmp:
            self.failUnlessRaises(ValueError, ps.write_graph,
                                  self.good_set, self.bad_edges, tmp)

    
class ReadPrimerTests(unittest.TestCase):

    def setUp(self):
        self.good_input = StringIO('''ATGC 1300 4000 3.1
GCCCCTA 1230 1300 1.1''')
        self.good_result = [ps.Primer(1, "ATGC", 4000, 1300),
                            ps.Primer(2, "GCCCCTA", 1300, 1230)]

        self.bad_input = StringIO('''ATGCGGC
TTTTCCC''')
        self.bad_handle = "primers.txt"
        
        
    def test_good_input(self):
        self.failUnlessEqual(ps.read_primers(self.good_input),
                             self.good_result)

    def test_bad_input(self):
        self.failUnlessRaises(ValueError, ps.read_primers,
                              self.bad_input)

    def test_bad_handle(self):
        self.failUnlessRaises(ValueError, ps.read_primers,
                              self.bad_handle)


class HeterodimerTests(unittest.TestCase):

    def setUp(self):
        self.primers = []
        self.primers.append(ps.Primer(1, "ATGCTC", 0, 0))
        # 4 contiguous
        self.primers.append(ps.Primer(2, "CAGCAT", 0, 0))
        # 3 contiguous, alt
        self.primers.append(ps.Primer(3, "GAGGTA", 0, 0))
        # 3 contiguous, alt 2
        self.primers.append(ps.Primer(4, "ATCGAC", 0, 0))
        # valid
        self.primers.append(ps.Primer(5, "TTCCAC", 0, 0))
                            

    def test_compatible_primers(self):
        edges = ps.test_pairs(self.primers[0:2], 3)
        self.assertEqual(edges, [])
        edges = ps.test_pairs([self.primers[0], self.primers[4]], 3)
        self.assertEqual(edges, [[1, 5]])


class FgBindLocationsTest(unittest.TestCase):

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
            primer = ps.Primer(1, "TATA", 0, 0)
            _, test_locs = ps.find_primer_locations(primer, filename)
        locs = sorted([0, 10, 1, 3, 5, 11, 2, 4, 6])
        self.assertEqual(test_locs, locs)


class FgBindDistanceTest(unittest.TestCase):

    def test_fg_bind_dist(self):
        self.line = "5 5 1 2 3"
        self.locations = {1:{'seq':'',
                             'loc':[1, 2, 3]},
                          2:{'seq':'',
                             'loc':[4, 5, 6]},
                          3:{'seq':'',
                             'loc':[7, 8, 10]}}
        primer_set, primers, max_dist, stdev = ps.fg_bind_distances(self.line, self.locations, ps.stdev)
        self.assertEqual(max_dist, 2)
        self.assertEqual(round(stdev, 4), 2.9345)
        self.assertEqual(primer_set, [1, 2, 3])
        



        
# def test_test_pairs():
#     pass

# def test_max_consecutive_binding():
#     pass
    
if __name__ == "__main__":
    unittest.main()
