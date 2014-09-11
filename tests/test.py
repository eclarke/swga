import tempfile
import unittest
from collections import namedtuple, OrderedDict
from StringIO import StringIO
import os
import swga
from swga.commands.filter_primers import filter_primers
from swga.commands.find_sets import find_sets







class HeterodimerTests(unittest.TestCase):
    '''
    Testing heterodimer checks
    '''
    def setUp(self):
        self.primers = []
        self.primers.append(swga.Primer(1, "ATGCTC", 0, 0, 0))
        # 4 contiguous
        self.primers.append(swga.Primer(2, "CAGCAT", 0, 0, 0))
        # 3 contiguous, alt
        self.primers.append(swga.Primer(3, "GAGGTA", 0, 0, 0))
        # 3 contiguous, alt 2
        self.primers.append(swga.Primer(4, "ATCGAC", 0, 0, 0))
        # valid
        self.primers.append(swga.Primer(5, "TTCCAC", 0, 0, 0))


    def test_compatible_primers(self):
        edges = swga.test_pairs(self.primers[0:2], 3)
        self.assertEqual(edges, [])
        edges = swga.test_pairs([self.primers[0], self.primers[4]], 3)
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
        locs = swga.find_locations(substr, test_str)
        self.assertEqual(locs, [0, 2, 4, 6])

    def test_find_primer_locations(self):
        test_str = ">TATATATAT>TATA"
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write(test_str)
            tmp.flush()
            filename = tmp.name
            primer = swga.Primer(1, "TATA", 0, 0, 0)
            _, test_locs = swga.find_primer_locations(primer, filename)
        locs = sorted([0, 10, 1, 3, 5, 11, 2, 4, 6])
        self.assertEqual(test_locs, locs)


if __name__ == "__main__":
    unittest.main()
