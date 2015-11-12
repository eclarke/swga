
# -*- coding: utf-8 -*-
"""
database.py

Defines the models and fields of the primer database using Peewee ORM. Also
contains a litany of helper functions for adding and retrieving stored data.

"""
import os
import swga
import json
from swga.utils import chunk_iterator
from swga.locate import revcomp
import swga.utils
import swga.score
import peewee as pw
try:
    from playhouse.shortcuts import ManyToManyField
except ImportError:
    from playhouse.fields import ManyToManyField

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
    found by using the PrimerSet intermediate table.
    '''
    _id = pw.IntegerField(null=True)
    seq = pw.CharField(primary_key=True)
    fg_freq = pw.IntegerField(default=0)
    bg_freq = pw.IntegerField(default=0)
    ratio = pw.FloatField(default=0.0)
    tm = pw.FloatField(null=True)
    _locations = pw.TextField(null=True)
    gini = pw.FloatField(null=True)
    active = pw.BooleanField(default=False)

    @staticmethod
    def exported_fields():
        fields = [
            'seq',
            'fg_freq',
            'bg_freq',
            'ratio',
            'gini',
            'tm']
        return fields

    def __repr__(self):
        rep_str = "Primer {0}:{1} (fg_freq:{2}, bg_freq:{3}, ratio:{4})"
        return rep_str.format(
            self.id, self.seq, self.fg_freq, self.bg_freq, self.ratio)

    def locations(self):
        if self._locations:
            return json.loads(self._locations)
        else:
            swga.error("No locations stored for " + str(self))

    def update_tm(self):
        self.tm = swga.melting.Tm(self.seq)

    def _update_locations(self, genome_fp):
        self._locations = json.dumps(
            swga.locate.binding_sites(self.seq, genome_fp))

    def _update_gini(self, genome_fp):
        chr_ends = swga.locate.chromosome_ends(genome_fp)
        locs = swga.locate.linearize_binding_sites([self], chr_ends)
        dists = swga.score.seq_diff(locs)
        self.gini = swga.stats.gini(dists)


class Set(SwgaBase):

    '''
    The sets table contains each set's metadata. The actual primers that belong
    in the set are found using the PrimerSet intermediate table.
    '''
    _id = pw.PrimaryKeyField()
    _hash = pw.IntegerField(unique=True, null=True)
    primers = ManyToManyField(Primer, related_name='sets')
    score = pw.FloatField()
    set_size = pw.IntegerField(null=True)
    bg_dist_mean = pw.FloatField(null=True)
    fg_dist_std = pw.FloatField(null=True)
    fg_dist_gini = pw.FloatField(null=True)
    fg_max_dist = pw.IntegerField(null=True)
    fg_dist_mean = pw.FloatField(null=True)
    scoring_fn = pw.TextField(null=True)

    def __repr__(self):
        attr_str = "; ".join(
            "{}:{}".format(k, v) for k, v
            in self.__dict__['_data'].items())
        return "Set: " + attr_str

    @staticmethod
    def exported_fields():
        fields = [
            '_id',
            'score',
            'set_size',
            'bg_dist_mean',
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
    found. Otherwise, it exits with an error (SystemExit).
    '''
    if db_fname is None:
        swga.error("Primer db name cannot be `None`: corrupt preferences.cfg?")
    elif db_fname == ":memory:":
        swga.warn("Creating in-memory primer database; this may not work.")
    elif not os.path.isfile(db_fname) and not create_if_missing:
        # Exits here
        swga.error(
            "Primer db not found at '%s': specify different filename or "
            "re-run `swga count`" % db_fname, exception=False
        )
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


def add_primers(primers, chunksize=199, add_revcomp=True):
    if add_revcomp:
        def mkrevcomp(p):
            p2 = dict(**p)
            p2['seq'] = revcomp(p['seq'])
            return p2
        primers += [mkrevcomp(p) for p in primers]
    chunk_iterator(
        primers,
        fn=lambda c: Primer.insert_many(c).execute(),
        n=chunksize,
        label="Updating database: ")


@db.atomic()
def add_set(_id, primers, **kwargs):
    if not primers:
        swga.error("Invalid primers for set")
    if isinstance(primers, pw.SelectQuery):
        nprimers = primers.count()
    else:
        nprimers = len(primers)
    if nprimers == 0:
        swga.error("Cannot have an empty set")
    _hash = hash(frozenset([p.seq for p in primers]))
    if Set.select(pw.fn.Count(Set._id)).where(Set._hash == _hash).scalar() > 0:
        return None
    s = Set.create(_id=_id, _hash=_hash, **kwargs)
    s.primers.add(primers)
    return s


def get_primers_for_set(set_id):
    set = Set.get(Set._id == set_id)
    return [p.seq for p in set.primers]


def get_primers_for_ids(pids):
    return list(Primer.select().where(Primer._id << pids).execute())
