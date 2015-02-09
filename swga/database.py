# -*- coding: utf-8 -*-
"""
database.py

Defines the models and fields of the primer database using Peewee ORM. Also
contains a litany of helper functions for adding and retrieving stored data.

"""
import os
import peewee as pw
import swga

# The primer database must be initialized before use
# ex: `db.init(db_fname)`
db = pw.SqliteDatabase(None)


def init_db(db_fname, create_if_missing=False):
    '''
    Initializes the database at the file path specified. 
    If `create_if_missing` is True, it will create the database if it can't be
    found. Otherwise, it throws an error.
    '''
    if not db_fname:
        swga.swga_error("Primer database name unspecified.")
    if create_if_missing and not os.path.isfile(db_fname):
        db.init(db_fname)
    elif not os.path.isfile(db_fname):
        swga.swga_error("Primer db not found at %s: specify different path or "
                        "re-run `swga count`" % db_fname)
    db.init(db_fname)
    return db


class SwgaBase(pw.Model):
    '''Specifies what database the inheriting models will use.'''
    class Meta:
        database = db
    

class Primer(SwgaBase):
    '''
    The primers table contains the sequence and metadata for each primer. Once
    set composition is determined, the sets that each primer belongs to can be
    found by using the Primer_Set intermediate table.
    '''
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

    def to_dict(self):
        return self.__dict__['_data']

    
class Set(SwgaBase):
    '''
    The sets table contains each set's metadata. The actual primers that belong
    in the set are found using the Primer_Set intermediate table.
    '''
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
        return "Set: "+"; ".join("{}:{}".format(k,v)
                                 for k,v in self.__dict__['_data_'])


class Primer_Set(SwgaBase):
    '''
    The "lookup table" enabling a many-to-many relationship between
    primers (which can belong to many sets) and sets (which contain many
    primers). 
    ''' 
    seq = pw.ForeignKeyField(Primer, related_name='sets', to_field='seq')
    set = pw.ForeignKeyField(Set, related_name='primers', to_field='sid')
    class Meta:
        indexes = (
            (('seq', 'set'), True),  # expects a tuple so needs a trailing ','
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


def update_in_chunks(itr, chunksize=100, model=Primer, show_progress=True,
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
        model.insert_many(p.to_dict() for p in chunk).upsert().execute()
    if isinstance(itr, pw.SelectQuery):
        itr = list(itr)    
    swga.core.chunk_iterator(itr, upsert_chunk, n=chunksize,
                             show_progress=show_progress,
                             label=label)
