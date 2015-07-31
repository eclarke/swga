# -*- coding: utf-8 -*-
'''
Handles bedgraph creation and export.
'''

import swga
import os
from types import *
from collections import Counter
import json


class BedGraph(object):
    """
    Creates a BedGraph-format file from a given set that illustrates the
    amount of coverage the primers in that set give over the foreground
    genome.
    """

    def __init__(self, set, fg_genome_fp, opts_str, window_size, step_size):
        '''
        args:
        set: the set of primers to export
        fg_genome_fp: the path to the foreground genome
        opts_str: the string to write after `track type=bedGraph`
        window_size: the size of the sliding window (<= size of each chromosome or record)
        step_size: the step size of the sliding window (must be smaller than window_size)
        '''
        super(BedGraph, self).__init__()
        assert isinstance(set, swga.database.Set)

        self.set = set
        self.fg_genome_fp = fg_genome_fp
        self.opts_str = opts_str
        self.window_size = int(window_size)
        self.step_size = int(step_size) if step_size else self.window_size / 5

    def _hits_per_record(self):
        '''
        Calculates the number of bases covered by primers in each record, and
        yields the record name, midpoint of the record, and number of hits as
        a tuple.
        '''
        record_ends = swga.locate.chromosome_ends(self.fg_genome_fp)
        for record_name, ends in record_ends.iteritems():
            record_length = ends[1] + 1

            # Check window size <= record_length and fix if not
            this_window_size = self.window_size
            if this_window_size > record_length:
                this_window_size = record_length
                swga.warn(
                    "In [{}]: window size larger than record; set to {}"
                    .format(record_name, this_window_size))

            # Check step size is compatible with window size and fix if not
            this_step_size = self.step_size
            if this_step_size > this_window_size:
                this_step_size = this_window_size
                swga.warn(
                    "In [{}]: step size larger than window size ({}), set to {}"
                    .format(record_name, this_window_size, this_step_size))

            # Count the number of primers that bind to any given nucleotide
            # in the current record
            counter = Counter()
            for primer in self.set.primers:
                k = len(primer.seq)
                locations = primer.locations()
                for l in locations[record_name]:
                    counter.update(Counter(xrange(l, l + k)))

            starting_positions = xrange(
                0, int(record_length - this_window_size), this_step_size)

            for start in starting_positions:
                end = start + this_window_size
                midpoint = (end + start) / 2
                # Add each base's count to get the number of bases covered
                hits = sum([counter[i] for i in xrange(start, end)])
                yield record_name, midpoint, hits

    def write(self, output_fp):
        """Writes the bedgraph to a file in a directory named after the set."""
        # Create the output folder if it doesn't exist already
        output_folder = swga.export._mk_folder(
            output_fp,
            self.fg_genome_fp,
            "set_%s" % self.set._id)

        bedgraph_fp = os.path.join(
            output_folder,
            "set_{}.bedgraph".format(self.set._id))

        with open(bedgraph_fp, 'wb') as bedgraph_file:
            typestr = "track type=bedGraph {}\n".format(self.opts_str)
            bedgraph_file.write(typestr)
            for record_name, midpoint, hits in self._hits_per_record():
                linestr = "{} {} {} {}\n".format(
                    record_name, midpoint, midpoint, hits)
                bedgraph_file.write(linestr)

        swga.message("Bedfile written to {}".format(bedgraph_fp))
