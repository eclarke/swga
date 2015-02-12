# -*- coding: utf-8 -*-
"""
database.py

Defines the models and fields of the primer database using Peewee ORM. Also
contains a litany of helper functions for adding and retrieving stored data.

"""
import os
import csv
import swga
import peewee as pw
from playhouse.shortcuts import ManyToManyField

# The primer database must be initialized before use
# ex: `db.init(db_fname)`
db = pw.SqliteDatabase(None)

class SwgaBase(pw.Model):
    '''Specifies what database the inheriting models will use.'''

    @classmethod
    def fieldnames(cls):
        '''Returns the fields defined in the model as a list of strings.'''
        return [name for name, attr in cls.__dict__.items() 
                if isinstance(attr, pw.FieldDescriptor)] 

    @classmethod
    def fields(cls):
        return dict((name, attr.field) for name, attr in cls.__dict__.items()
                    if isinstance(attr, pw.FieldDescriptor))

    def to_dict(self):
        return self.__dict__['_data']
    
    class Meta:
        database = db
    

class Primer(SwgaBase):
    '''
    The primers table contains the sequence and metadata for each primer. Once
    set composition is determined, the sets that each primer belongs to can be
    found by using the Primer_Set intermediate table.
    '''
    _id = pw.IntegerField(null=True)
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

    @staticmethod
    def exported_fields():
        fields = [
            'seq',
            'fg_freq',
            'bg_freq',
            'ratio',
            'tm'
        ]
        return fields

    
    
class Set(SwgaBase):
    '''
    The sets table contains each set's metadata. The actual primers that belong
    in the set are found using the PrimerSet intermediate table.
    '''
    _id = pw.PrimaryKeyField()
    primers = ManyToManyField(Primer, related_name='sets')
    score = pw.FloatField()
    set_size = pw.IntegerField(null=True)
    bg_ratio = pw.FloatField(null=True)
    fg_dist_std = pw.FloatField(null=True)
    fg_dist_gini = pw.FloatField(null=True)
    fg_max_dist = pw.IntegerField(null=True)
    fg_dist_mean = pw.FloatField(null=True)
    scoring_fn = pw.TextField(null=True)

    def __repr__(self):
        return "Set: "+"; ".join("{}:{}".format(k,v)
                                 for k,v in self.__dict__['_data'].items())

    @staticmethod
    def exported_fields():
        fields = [
            'score',
            'set_size',
            'bg_ratio',
            'fg_max_dist',
            'fg_dist_mean',
            'fg_dist_std',
            'fg_dist_gini',
            'scoring_fn',
            'primers'
        ]
        return fields
            

PrimerSet = Set.primers.get_through_model()
    

def init_db(db_fname, create_if_missing=False):
    '''
    Initializes the database at the file path specified. 
    If `create_if_missing` is True, it will create the database if it can't be
    found. Otherwise, it throws an error.
    '''
    if not db_fname:
        swga.swga_error("Primer database name unspecified.")
    if db_fname == ":memory:":
        db.init(db_fname)
    elif create_if_missing and not os.path.isfile(db_fname):
        db.init(db_fname)
    elif not os.path.isfile(db_fname):
        swga.swga_error("Primer db not found at %s: specify different path or "
                        "re-run `swga count`" % db_fname)
    db.init(db_fname)
    return db


def create_tables(drop=True):
    '''
    Creates the tables in the database. If `drop`=True, attempts to safely drop
    the tables before creating them.
    '''
    tbl_list = [Primer, Set, PrimerSet]
    if drop:
        db.drop_tables(tbl_list, safe=True)  
    db.create_tables(tbl_list, safe=True)


def add_set(_id, primers, **kwargs):
    if not primers:
        swga.swga_error("Invalid primers for set")
    if isinstance(primers, pw.SelectQuery):
        nprimers = primers.count()
    else:
        nprimers = len(primers)
    if nprimers == 0:
        swga.swga_error("Cannot have an empty set")
    Set.delete().where(Set._id == _id)
    s = Set.create(_id=_id, **kwargs)
    s.primers.add(primers)
    return s


def get_primers_for_set(set_id):
    set = Set.get(Set._id == set_id)
    return [p.seq for p in set.primers]


def get_primers_for_ids(pids):
    return list(Primer.select().where(Primer._id << pids).execute())


def update_in_chunks(itr, chunksize=100, show_progress=True,
                     label=None):
    '''
    Inserts or updates records in database in chunks of a given size.

    Arguments:
    - itr: a list or other iterable containing records in the primer db that
           have a to_dict() method
    - chunksize: the size of the chunk. Usually has to be 999/(number of fields)
    - model: the table in the db to update
    - show_progress, label: passed to progress.bar
    '''
    def upsert_chunk(chunk):
        seqs = [p.seq for p in chunk]
        Primer.delete().where(Primer.seq << seqs).execute()
        Primer.insert_many(p.to_dict() for p in chunk).execute()
    if isinstance(itr, pw.SelectQuery):
        itr = list(itr)
    swga.core.chunk_iterator(itr, upsert_chunk, n=chunksize,
                             show_progress=show_progress,
                             label=label)

    

