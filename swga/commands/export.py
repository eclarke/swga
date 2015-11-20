# -*- coding: utf-8 -*-
'''
Exports primers or sets from the database in a variety of formats.

Examples:

Get set metadata and primers for top 10 sets by score
$ swga export sets --order_by score --descending --limit 10

Get metadata for these four sets
$ swga export sets --ids 1 4 6 12

Get 50 primers with the lowest bg binding frequency
$ swga export primers --order_by bg_freq --limit 50

Get these primers' metadata
$ swga export primers --seq ATGGGCCA GCCATT

Get a set as a BED file for IGV
$ swga export bed_file --set_id 5

Get the Lorenz curve for a set (a visualization of the Gini index)
$ swga export lorenz --set_id 5
'''
import collections
import csv
import os

import swga.database as database
import swga.locate
import swga.stats
from swga.commands._command import Command
from swga.database import Primer, Set
from swga.export import BedGraph, BedFile


class Export(Command):

    def run(self):
        self.header = not self.no_header
        if self.limit < 0:
            self.limit = None
        if self.what in ['set', 'sets']:
            sets = self.get_items(Set)
            export(Set, sets, self.output, self.header)

        if self.what in ['primer', 'primers']:
            primers = self.get_items(Primer)
            export(Primer, primers, self.output, self.header)

        if self.what == 'lorenz':
            sets = self.get_items(Set)
            export_lorenz(sets, self.output, self.fg_genome_fp, self.header)

        if "bed" in self.what:
            outpath = self.output_folder if self.output_folder else os.getcwd()
            sets = self.get_items(Set)
            if self.what == "bedfile":
                for set in sets:
                    swga.message(
                        "Exporting set {} and its primers as bedfiles..."
                        .format(set._id))
                    bedfile = BedFile(set, self.fg_genome_fp)
                    bedfile.write(outpath)
            elif self.what == "bedgraph":
                for set in sets:
                    swga.message("Exporting set {} as bedgraph..."
                                 .format(set._id))
                    bedgraph = BedGraph(
                        set=set,
                        fg_genome_fp=self.fg_genome_fp,
                        opts_str=self.opts_str,
                        window_size=self.window_size,
                        step_size=self.step_size)
                    bedgraph.write(outpath)

    def get_items(self, model):
        '''
        Retrieves the rows from the model given by the arguments (basically an
        adaptor around very simple SQL queries.)
        '''
        ids = self.ids
        order_by = self.order_by
        limit = self.limit
        descending = self.descending

        if ids and (order_by or limit):
            swga.warn("ID(s) specified: ignoring --order_by/--limit options.")
            order_by = limit = None

        validate_order_field(order_by, model)

        if ids:
            targets = model.select().where(model._id << ids)
            for target in targets:
                yield target
        else:
            query = model.select()
            if order_by:
                if descending:
                    query = query.order_by(model.fields()[order_by].desc())
                else:
                    query = query.order_by(model.fields()[order_by])
            if limit:
                query = query.limit(limit)
            targets = query
            for target in targets:
                yield target


def export(model, rows, outfile, header=True):
    '''
    Writes rows from a model (i.e. table) in a tab-delimited format to a file
    specified by outfile. Relies on the custom `exported_fields` function of the
    Primer and Set models.
    '''
    exported_fields = model.exported_fields()

    writer = csv.DictWriter(
        outfile,
        exported_fields,
        delimiter="\t",
        extrasaction='ignore',
        lineterminator='\n'
    )

    if header:
        writer.writeheader()

    for row in rows:
        d = {}
        for field_name in row.exported_fields():
            field = getattr(row, field_name)
            if (not isinstance(field, basestring) and
                    isinstance(field, collections.Iterable)):
                try:
                    field = ",".join([p.seq for p in field])
                except AttributeError:
                    field = ",".join([s._id for s in field])
            d[field_name] = field
        writer.writerow(d)


def export_lorenz(sets, outfile, fg_genome_fp, header=True):
    '''Exports the empirical Lorenz curve used for the given sets.'''
    writer = csv.DictWriter(
        outfile,
        fieldnames=["SetID", "CDF"],
        delimiter="\t",
        lineterminator='\n'
    )

    if header:
        writer.writeheader()

    for set in sets:
        # Get the distances between each primer binding site
        primer_seqs = database.get_primers_for_set(set._id)
        primers = list(Primer.select().where(Primer.seq << primer_seqs).execute())
        chr_ends = swga.locate.chromosome_ends(fg_genome_fp)
        binding_sites = swga.locate.linearize_binding_sites(primers, chr_ends)
        distances = swga.stats.seq_diff(binding_sites)

        lorenz = swga.stats.lorenz(distances)
        lorenz_str = ",".join(str(d) for d in lorenz)
        writer.writerow({'SetID': set._id, 'CDF': lorenz_str})


def validate_order_field(field, model):
    '''Ensures the given field exists in the model.'''
    if field and field not in model.fields():
        swga.error(
            "Cannot order by '{}'. Valid choices are {}"
            .format(field, ", ".join(Primer.fields())))
