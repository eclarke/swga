# -*- coding: utf-8 -*-
'''
Handles bedfile creation and export.
'''

import swga
import os


class BedFile(object):
    """Creates and writes BedFiles from a set of primers."""
    def __init__(self, set, fg_genome_fp):
        super(BedFile, self).__init__()
        assert isinstance(set, swga.database.Set)

        self.set = set
        self.fg_genome_fp = fg_genome_fp

    def _primer_sites(self, primer):
        seq = primer.seq
        for record_name, locations in primer.locations_dict().iteritems():
            for location in locations:
                yield record_name, location, location + len(seq)

    def write(self, output_fp):
        '''
        Writes n+1 BED-formatted files to a folder, where n is the number of
        primers in a set. One file is written for each primer, and then one
        file is written that combines all the primers together into one set.
        '''
        # Create the output folder if it doesn't exist already
        output_folder = swga.export._mk_folder(
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
