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

import swga.database as database
from swga import message, warn, swga_error
from swga.database import Primer, Set
from swga.commands import Command

def main(argv, cfg_file):
    cmd = Command('export', cfg_file=cfg_file)
    cmd.parse_args(argv)
    
    what = cmd.what
    database.init_db(cmd.primer_db)
    if what in ['set', 'sets']:
        for _set in export(Set, cmd.ids, cmd.order_by, cmd.limit, cmd.descending):
            cmd.output.write(_set+'\n')
    if what in ['primer', 'primers']:
        export(Primer, cmd.ids, cmd.order_by, cmd.limit)
    if what == "bedfile":
        warn("Not yet implemented.")
        pass


def validate_args(export_fn):
    def validate(model, ids, order_by, limit, descending):
        if ids and (order_by or limit):
            warn("List of IDs specified, ignoring order_by and limit parameters")
            order_by = limit = None
        validate_order_field(order_by, model)
        return export_fn(model, ids, order_by, limit, descending)
    return validate


@validate_args
def export(model, ids=None, order_by=None, limit=None, descending=False):
    if ids:
        print ids
        targets = model.select().where(model._id << ids)
        print targets
    else:
        query = model.select()
        if order_by:
            if descending:
                print descending
                query = query.order_by(model.fields()[order_by].desc())
                print query
            else:
                query = query.order_by(model.fields()[order_by])
        if limit:
            query = query.limit(limit)
        targets = query
        for target in targets:
            yield str(target)
    


def validate_order_field(field, model):
    if field and field not in model.fields():
        swga_error("Cannot order by '{}'. Valid choices are {}"
                   .format(field, ", ".join(Primer.fields())))

