# -*- coding: utf-8 -*-
"""Classes that define a swga workspace.

- SwgaWorkspace:
    A database connection where primers and sets are stored.
- SwgaModel:
    The base class for records (primers, sets) that are stored in the database.
"""

import json
import contextlib

import semantic_version as semver
import peewee as pw
try:
    from playhouse.shortcuts import ManyToManyField
except ImportError:
    from playhouse.fields import ManyToManyField
from playhouse.sqlite_ext import SqliteExtDatabase

from swga import error, meta, __version__
import locate
import melting
import stats


class SwgaWorkspace(SqliteExtDatabase):

    """Extends a SqliteExtDatabase to add workspace metadata getting/setting."""

    @property
    def metadata(self):
        m = _metadata.get()
        m = {key: getattr(m, key) for key in m.fieldnames() if key != 'id'}
        _m = meta(**m)
        return _m

    @metadata.setter
    def metadata(self, values):
        _metadata.delete().execute()
        _metadata.insert(db_name=self.database, **values).execute()

    def create_tables(self, safe=True):
        """Create the tables that subclass SwgaModel."""
        super(SwgaWorkspace, self).create_tables(_tables, safe=safe)

    def reset_sets(self):
        self.drop_tables([Set, PrimerSet])
        super(SwgaWorkspace, self).create_tables([Set, PrimerSet])

    def reset_primers(self):
        self.drop_tables([Set, PrimerSet, Primer])
        super(SwgaWorkspace, self).create_tables([Set, PrimerSet, Primer])

    def check_version(self, version):
        """Check the version of the database and compare it to the swga version.

        If the two versions are incompatible, raise a SystemExit.
        """
        try:
            db_ver = self.metadata.version
            db_ver = semver.Version(self.metadata.version)
        except pw.OperationalError:
            db_ver = "<NA>"
        ver = semver.Version(version)
        spec = semver.Spec('=={}.{}'.format(ver.major, ver.minor))
        if db_ver not in spec:
            error(
                "This workspace was created with a different version of swga.\n"
                "  Workspace version: {}\n"
                "  swga version:      {}\n"
                "Please re-initialize the workspace with `swga init` or use a "
                "different version of swga."
                .format(db_ver, ver),
                exception=False,
                wrap=False
            )

_db = SwgaWorkspace(None)


class SwgaModel(pw.Model):

    """Specifies database location and model-generalized methods."""

    @classmethod
    def fieldnames(cls):
        """Return the fields defined in the model as a list of strings."""
        return [name for name, attr in cls.__dict__.items()
                if isinstance(attr, pw.FieldDescriptor)]

    @classmethod
    def fields(cls):
        """Return a dictionary of field/values."""
        return dict((name, attr.field) for name, attr in cls.__dict__.items()
                    if isinstance(attr, pw.FieldDescriptor))

    def to_dict(self):
        """Return internal data of the model as a dictionary."""
        return self.__dict__['_data']

    class Meta:
        database = _db


class Primer(SwgaModel):

    """
    The primers table contains the sequence and metadata for each primer.

    Once set composition is determined, the sets that each primer belongs to can
    be found by using the PrimerSet intermediate table.
    """

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
        rep_str = "Primer {0}:{1}".format(self._id, self.seq)
        return rep_str

    @property
    def locations(self):
        if self._locations:
            return json.loads(self._locations)
        else:
            error("No locations stored for " + str(self))

    def update_tm(self):
        self.tm = melting.temp(self.seq)

    def _update_locations(self, genome_fp):
        self._locations = json.dumps(
            locate.binding_sites(self.seq, genome_fp))

    def _update_gini(self, genome_fp):
        chr_ends = locate.chromosome_ends(genome_fp)
        locs = locate.linearize_binding_sites([self], chr_ends)
        dists = stats.seq_diff(locs)
        self.gini = stats.gini(dists)


class Set(SwgaModel):

    """
    The sets table contains each set's metadata. The actual primers that belong
    in the set are found using the PrimerSet intermediate table.
    """

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

    def primer_seqs(self):
        return [p.seq for p in self.primers]

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

    @staticmethod
    @_db.atomic()
    def add(_id, primers, **kwargs):
        """Add a set with the specified primers to the database.

        If the same set exists in the database already, this does not add it again.
        :param _id: the set ID. For user-defined sets, this should be negative.
        :param primers: a Primers object defining the set to add
        :param kwargs: additional arguments to Set.create
        :return: the new Set object and True if the set was added or False if it was not
        """
        try:
            primers = primers.primers
        except AttributeError:
            pass

        if ((not primers) or
            (isinstance(primers, pw.SelectQuery) and primers.count == 0) or
            (len(primers) == 0)):
            raise ValueError("Cannot add an empty set.")

        # We convert the primers into a hashable set to check and see if a set
        # with the same primers already exists in the db
        _hash = hash(frozenset([p.seq for p in primers]))
        s, created = Set.get_or_create(_hash=_hash, defaults=dict(_id=_id, **kwargs))
        if created:
            s.primers.add(primers)
        return s, created

PrimerSet = Set.primers.get_through_model()


class _metadata(SwgaModel):

    """Workspace metadata (internal)."""

    db_name = pw.TextField()
    version = pw.TextField(null=True)
    fg_file = pw.TextField(null=True)
    bg_file = pw.TextField(null=True)
    ex_file = pw.TextField(null=True)
    fg_length = pw.IntegerField(null=True)
    bg_length = pw.IntegerField(null=True)


# Uses new-class introspection to find all the models that subclass SwgaModel,
# then adds them to the _tables list for use in other functions. This is module-
# level because we need it to be available to classes within the module.
_tables = [cls for cls in vars()['SwgaModel'].__subclasses__()]
_tables.append(PrimerSet)


@contextlib.contextmanager
def connection(db_name):
    """Return a connection context to the database.

    Ensures the database is initialized before use, and closed when finished.
    """

    _db.init(db_name)
    _db.connect()
    yield _db
    _db.close()

