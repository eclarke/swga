import os
from types import *
from collections import Counter
from swga import (warn, message)
import swga.database as database
import swga.locate as locate


def _mk_folder(folder, fg_genome_fp, setstr):
    fg_name = ".".join(
        os.path.basename(fg_genome_fp).split(".")[:-1]) + "_export"
    output_folder = os.path.join(folder, fg_name, setstr)
    if not os.path.isdir(output_folder):
            os.makedirs(output_folder)
    return output_folder


class BedFile(object):
    """Creates and writes BedFiles from a set of primers."""
    def __init__(self, set, fg_genome_fp):
        super(BedFile, self).__init__()
        assert isinstance(set, database.Set)

        self.set = set
        self.fg_genome_fp = fg_genome_fp

    def _primer_sites(self, primer):
        seq = primer.seq
        for record_name, locations in primer.locations.iteritems():
            for location in locations:
                yield record_name, location, location + len(seq)

    def write(self, output_fp):
        '''
        Writes n+1 BED-formatted files to a folder, where n is the number of
        primers in a set. One file is written for each primer, and then one
        file is written that combines all the primers together into one set.
        '''
        # Create the output folder if it doesn't exist already
        output_folder = _mk_folder(
            output_fp,
            self.fg_genome_fp,
            "set_%s" % self.set._id)

        whole_set_fp = os.path.join(output_folder, 'whole_set.bed')

        with open(whole_set_fp, 'wb') as whole_set_file:
            whole_set_file.write("track name=Set_{}\n".format(self.set._id))
            for primer in self.set.primers:
                seq = primer.seq
                primer_fp = os.path.join(output_folder, seq + ".bed")
                with open(primer_fp, 'wb') as primer_file:
                    primer_file.write("track name={}\n".format(seq))
                    for record_name, start, end in self._primer_sites(primer):
                        recordstr = "{} {} {}\n".format(
                            record_name, start, end)
                        primer_file.write(recordstr)
                        whole_set_file.write(recordstr)
        message("Bedfiles written to "+output_folder)


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
        assert isinstance(set, database.Set)

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
        record_ends = locate.chromosome_ends(self.fg_genome_fp)
        for record_name, ends in record_ends.iteritems():
            record_length = ends[1] + 1

            # Check window size <= record_length and fix if not
            this_window_size = self.window_size
            if this_window_size > record_length:
                this_window_size = record_length
                warn(
                    "In [{}]: window size larger than record; set to {}"
                    .format(record_name, this_window_size))

            # Check step size is compatible with window size and fix if not
            this_step_size = self.step_size
            if this_step_size > this_window_size:
                this_step_size = this_window_size
                warn(
                    "In [{}]: step size larger than window size ({}), set to {}"
                    .format(record_name, this_window_size, this_step_size))

            # Count the number of primers that bind in the current record
            counter = Counter()
            for primer in self.set.primers:
                counter.update(Counter(primer.locations[record_name]))

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
        output_folder = _mk_folder(
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

        message("Bedgraph written to {}".format(bedgraph_fp))