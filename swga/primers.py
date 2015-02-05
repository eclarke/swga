# -*- coding: utf-8 -*-
"""primers.py

Functions for parsing primers from a file.

"""
from __future__ import with_statement, division
import subprocess
import csv
import os
import re
import struct
import swga
import swga.resources as resources
import peewee as pw
from pyfaidx import Fasta


# The primer database must be initialized before use
# ex: `db.init(db_fname)`
db = pw.SqliteDatabase(None)


def init_db(db_fname, create_if_missing=False):
    fp = db_fname
    if not os.path.isfile(fp) and create_if_missing:
        db.init(fp)
    elif not os.path.isfile(fp):
        swga.swga_error("Primer db not found at %s: specify different path or "
                        "re-run `swga count`" % fp)
    db.init(fp)
    return db


class SwgaBase(pw.Model):
    class Meta:
        database = db
    

class Primer(SwgaBase):
    pid = pw.IntegerField(null=True)
    seq = pw.CharField(primary_key=True)
    fg_freq = pw.IntegerField(default=0)
    bg_freq = pw.IntegerField(default=0)
    ratio = pw.FloatField(default=0.0)
    tm = pw.FloatField(null=True)
    locations = pw.TextField(null=True)
    active = pw.BooleanField(default=False)

    def __repr__(self):
        rep_str = "Primer {0}:{1} (fg_freq:{2}, bg_freq:{3}, ratio:{4})"
        return rep_str.format(
            self.id, self.seq, self.fg_freq, self.bg_freq, self.ratio)


class Set(SwgaBase):
    sid = pw.PrimaryKeyField()
    pids= pw.TextField(unique=True)
    score = pw.FloatField()
    set_size = pw.IntegerField(null=True)
    bg_ratio = pw.FloatField(null=True)
    fg_dist_std = pw.FloatField(null=True)
    fg_dist_gini = pw.FloatField(null=True)
    fg_max_dist = pw.IntegerField(null=True)
    fg_dist_mean = pw.FloatField(null=True)
    scoring_fn = pw.TextField(null=True)

    def __repr__(self):
        return "Set: "+"; ".join("{}:{}".format(k,v) for k,v in self.__dict__)


class Primer_Set(SwgaBase):
    seq = pw.ForeignKeyField(Primer, related_name='sets', to_field='seq')
    set = pw.ForeignKeyField(Set, related_name='primers', to_field='sid')
    class Meta:
        indexes = (
            (('seq', 'set'), True),
        )


def create_tables(drop=True):
    if drop:
        db.drop_tables([Primer, Set, Primer_Set], safe=True)        
    db.create_tables([Primer, Set, Primer_Set], safe=True)


def add_set(primers, **kwargs):
    try:
        s = Set.create(**kwargs)
        for primer in primers:
            Primer_Set.create(seq=primer.seq, set=s)
        return set
    except pw.IntegrityError:
        pass


def get_primers_for_set(set_id):
    set = Set.get(Set.sid == set_id)
    return [p.seq.seq for p in set.primers]


def get_primers_for_ids(pids):
    return list(Primer.select().where(Primer.pid << pids).execute())


def count_kmers(k, genome_fp, cwd, threshold=1):
    assert isinstance(threshold, int)
    dsk = resources.get_dsk()
    genome = genome_fp.split(os.sep).pop()
    out = '%s-%dmers' % (genome, k)
    outfile = os.path.join(cwd, out + '.solid_kmers_binary')
    if os.path.isfile(outfile):
        swga.message("Binary kmer file already found at %s, skipping..."
                     % outfile)
    else:
        cmdstr = ("{dsk} {genome_fp} {k} -o {out} -t {threshold}"
                  .format(**locals()))
        swga.message("In {cwd}:\n> {cmdstr}".format(**locals()))
        try:
            subprocess.check_call(cmdstr, shell=True, cwd=cwd)
        except:
            if os.path.isfile(outfile):
                os.remove(outfile)
            raise
    primers = dict((kmer, freq)
                   for kmer, freq in parse_kmer_binary(outfile)
                   if freq >= threshold)
    return primers
        

def parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        try:
            kmer_nbits = struct.unpack('i', f.read(4))[0]
            k = struct.unpack('i', f.read(4))[0]
        except struct.error:
            if os.path.isfile(fp):
                os.remove(fp)
            raise
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits // 8),
                                            f.read(kmer_nbits // 8))
                freq = struct.unpack('I', f.read(4))[0]

                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[(kmer_binary[i//4] >> (2 * (i%4))) % 4] + kmer

                yield kmer, freq

        except struct.error:
            pass


def get_primer_locations(seq, genome_fp):
    genome = Fasta(genome_fp)
    len_so_far = 0
    locations = []
    for chromosome in genome:
        chr_len = len(chromosome)
        fwd_locs = find_locations(str(seq), str(chromosome))
        rev_locs = find_locations(str(seq)[::-1], str(chromosome))
        locations += [l + len_so_far for l in fwd_locs + rev_locs]
        len_so_far += chr_len
    if not locations:
        print str(chromosome[1:100])
        print str(seq)
        raise ValueError("No binding sites!")
    return locations


def get_chromosome_ends(genome_fp):
    genome = Fasta(genome_fp)
    len_so_far = 0
    chr_ends = []
    for chromosome in genome:
        chr_len = len(chromosome)
        chr_ends += [len_so_far, chr_len + len_so_far]
    return chr_ends


def find_locations(substring, string):
    '''
    Very fast way of finding overlapping substring locations in a
    (potentially large) string.
    '''
    locations = []
    start = 0
    # Assumes string is all upper-case
    substring = substring.upper()
    while True:
        start = string.find(substring, start) + 1
        if start > 0:
            locations.append(start-1)
        else:
            return locations


def read_primer_file(infile, echo_input=False):
    '''
    Calls parse_primer() on each line of the input file. If a malformed line is
    found, will skip parsing and (optionally) output a warning message.

    Arguments:
    infile: an open file handle to read from

    Returns: A list of Primer objects
    '''
    rows = csv.DictReader((f for f in infile if not f.startswith('#')),
                          delimiter='\t',
                          fieldnames=['seq', 'fg_freq', 'bg_freq', 'ratio'], 
                          restkey='extra')
    primers = [Primer(id=i+1, **row) for i, row in enumerate(rows)]
    return primers


def read_primer_list(plist):
    '''
    Reads in a list of primers, one per line, and returns the corresponding
    records from the primer database.
    '''
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in plist]
    return list(Primer.select().where(Primer.seq << seqs).execute())


def write_primer_file(primers, out_handle, header=True):
    '''Writes each primer in primers to a line.

    Arguments:
    - primers: a list of primers
    - out_handle: an open file handle to write to
    - header: write a header?
    '''
    fields = ["seq", "fg_freq", "bg_freq", "ratio"]
    writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter='\t', 
                            lineterminator='\n',
                            extrasaction='ignore')
    if header:
        out_handle.write("#")
        writer.writeheader()
    [writer.writerow(p.__dict__) for p in primers]


