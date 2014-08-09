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
        self.failUnlessRaises(AttributeError, ps.read_primers,
                              self.bad_handle)


class HeterodimerTests(unittest.TestCase):

    def setUp(self):
        self.primer1 = ps.Primer(1, "ATGCTC", 0, 0)
        # 3 contiguous
        self.primer2 = ps.Primer(2, "TACAAC", 0, 0)
        # 3 contiguous, alt
        self.primer3 = ps.Primer(3, "ATGGAG", 0, 0)
        # 3 contiguous, alt 2
        self.primer4 = ps.Primer(4, "ATCGAC", 0, 0)
        # valid
        self.primer5 = ps.Primer(5, "TTCCAC", 0, 0)
                            

    def test_compatible_primers(self):
        

        
# def test_test_pairs():
#     pass

# def test_max_consecutive_binding():
#     pass
    
if __name__ == "__main__":
    unittest.main()
