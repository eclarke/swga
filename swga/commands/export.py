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
import json
import os
from collections import defaultdict, Counter

import swga.database as database
import swga.locate
import swga.stats
from swga import warn, swga_error
from swga.commands import Command
from swga.database import Primer, Set
from swga.clint.textui import progress

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
        primers = get_items(Primer, cmd.ids, cmd.order_by, cmd.limit, cmd.descending)
        export(Primer, primers, cmd.output, header)

    if "bed" in what:
        outpath = cmd.output_folder if cmd.output_folder else os.getcwd()
        sets = get_items(Set, cmd.ids, cmd.order_by, cmd.limit, cmd.descending)        
        if what == "bedfile":
            for set in sets:
                swga.message(
                    "Exporting set {} and associated primers as bedfiles..."
                    .format(set._id))
                export_bedfiles(set, cmd.fg_genome_fp, outpath)
        elif what == "bedgraph":
            for set in sets:
                swga.message("Exporting set {} as bedgraph...".format(set._id))
                export_bedgraph(
                    set=set,
                    fg_genome_fp=cmd.fg_genome_fp,
                    outpath=outpath,
                    opts_str=cmd.opts_str,
                    window_size=cmd.window_size,
                    step_size=cmd.step_size)
            


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


def _mk_folder_name(folder, fg_genome_fp, setstr):
    fg_name = ".".join(os.path.basename(fg_genome_fp).split(".")[:-1]) + "_export"
    return os.path.join(folder, fg_name, setstr)
                               
        
def export_bedfiles(set, fg_genome_fp, outpath):
    '''
    Writes n+1 BED-formatted files to a folder, where n is the number of primers
    in a set. One file is written for each primer, and then one file is written
    that combines all the primers together into one set.
    '''

    output_folder = _mk_folder_name(outpath, fg_genome_fp, "set_%s" % set._id)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
        
    whole_set_fp = os.path.join(output_folder, 'whole_set.bed')

    with open(whole_set_fp, 'wb') as whole_set_file:
        whole_set_file.write("track name=Set_{}\n".format(set._id))
        # Each location of each primer gets its own line in this file 
        for primer in set.primers:
            seq = primer.seq
            primer_fp = os.path.join(output_folder, seq)
            with open(primer_fp, 'wb') as primer_file:
                primer_file.write("track name={}\n".format(seq))
                for record in json.loads(primer.locations):
                    record_name = record.split('|')[0].strip()
                    for location in primer.locations[record]:
                        record_string = "{} {} {}\n".format(
                            record_name, location, location+len(primer.seq))
                        primer_file.write(record_string)
                        whole_set_file.write(record_string)


def export_bedgraph(set, fg_genome_fp, outpath, opts_str, window_size, step_size):

    output_folder = _mk_folder_name(outpath, fg_genome_fp, "set_%s" % set._id)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
        
    step_size = step_size if step_size else int(window_size/5)
    
    chr_ends = swga.locate.chromosome_ends(fg_genome_fp)
    records = chr_ends.keys()
    chr_lengths = dict((r, chr_ends[r][1] + 1) for r in records)
    bedgraph_fp = os.path.join(output_folder, "{}.bedgraph".format(set._id))
    with open(bedgraph_fp, 'wb') as bedgraph_file:
        typestr = "track type=bedGraph {}\n".format(opts_str)
        bedgraph_file.write(typestr)
        for record in records:
            record_name = record.split('|')[0].strip()
            chr_len = chr_lengths[record]

            this_window_size = window_size
            if window_size > chr_len:
                this_window_size = chr_len
                swga.warn("In [{}]: window size larger than record; set to {}"
                          .format(record_name, this_window_size))

            this_step_size = step_size
            if this_step_size > this_window_size:
                this_step_size = this_window_size / 5    
                swga.warn("Step size larger than window; set to {}"
                          .format(record_name, this_step_size))

            this_window_size = int(this_window_size)
            this_step_size = int(this_step_size)
            
            # Counter objects tally the unique items assigned to them.
            # This one tallies the number of primers that bind to any given
            # nucleotide in the FASTA record.
            # Adding counters together increments common items and merges
            # missing items.
            counter = Counter()
            for primer in set.primers:
                k = len(primer.seq)
                for l in json.loads(primer.locations)[record]:
                    counter.update(Counter(xrange(l, l + k)))

            window_iter = xrange(0, int(chr_len - this_window_size), this_step_size)
            for start in progress.bar(
                    window_iter,
                    expected_size=int(chr_len/this_step_size),
                    label=record_name):
                end = start + this_window_size
		midpoint = (end + start)/2
                # number of bases bound ('hit')
                nhits = sum([counter[i] for i in xrange(start, end)])
                linestr = "{} {} {} {}\n".format(
                    record_name, midpoint, midpoint, nhits)
                bedgraph_file.write(linestr)
            swga.message("Bedfile written to {}".format(bedgraph_fp))
            
                
            
