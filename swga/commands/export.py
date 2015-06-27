# -*- coding: utf-8 -*-
'''
Exports primers or sets from the database.
Basically a thin wrapper around a bunch of SQLite commands, and some export
formatting.

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
'''
import collections
import csv
import os

import swga.database as database
import swga.locate
import swga.stats
from swga import warn, swga_error
from swga.commands import Command
from swga.database import Primer, Set
from swga.export import BedGraph, BedFile


def main(argv, cfg_file):
    cmd = Command('export', cfg_file=cfg_file)
    cmd.parse_args(argv)
    header = not cmd.no_header
    what = cmd.what
    database.init_db(cmd.primer_db)

    if what in ['set', 'sets']:
        sets = get_items(Set, cmd.ids, cmd.order_by, cmd.limit, cmd.descending)
        export(Set, sets, cmd.output, header)

    if what in ['primer', 'primers']:
        primers = get_items(
            Primer, cmd.ids, cmd.order_by, cmd.limit, cmd.descending)
        export(Primer, primers, cmd.output, header)

    if "bed" in what:
        outpath = cmd.output_folder if cmd.output_folder else os.getcwd()
        sets = get_items(Set, cmd.ids, cmd.order_by, cmd.limit, cmd.descending)
        if what == "bedfile":
            for set in sets:
                swga.message(
                    "Exporting set {} and associated primers as bedfiles..."
                    .format(set._id))
                bedfile = BedFile(set, cmd.fg_genome_fp)
                bedfile.write(outpath)
        elif what == "bedgraph":
            for set in sets:
                swga.message("Exporting set {} as bedgraph...".format(set._id))
                bedgraph = BedGraph(
                    set=set,
                    fg_genome_fp=cmd.fg_genome_fp,
                    opts_str=cmd.opts_str,
                    window_size=cmd.window_size,
                    step_size=cmd.step_size)
                bedgraph.write(outpath)


def get_items(model, ids=None, order_by=None, limit=None, descending=False):
    '''
    Retrieves the rows from the model given by the arguments (basically an
    adaptor around very simple SQL queries.)
    '''

    if ids and (order_by or limit):
        warn("ID(s) specified: ignoring --order_by and --limit options.")
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
        extrasaction='ignore')

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


def validate_order_field(field, model):
    '''Ensures the given field exists in the model.'''
    if field and field not in model.fields():
        swga_error(
            "Cannot order by '{}'. Valid choices are {}"
            .format(field, ", ".join(Primer.fields())))
