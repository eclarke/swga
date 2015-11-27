
# -*- coding: utf-8 -*-
"""
database.py

Defines the models and fields of the primer database using Peewee ORM. Also
contains a litany of helper functions for adding and retrieving stored data.

"""
import json
import os
from contextlib import contextmanager
import click
from swga import (error, warn, meta, __version__)
from swga.utils import chunk_iterator
import swga.stats as stats
import swga.locate as locate
import melting
import semantic_version as semver


@contextmanager
def connection(db_name, create_if_missing=False):
    if db_name != ':memory:':
        if not os.path.isfile(db_name) and not create_if_missing:
            error(
                '{} does not appear to be a workspace (database not found). '
                'Run `swga init` to initialize this directory as a workspace.'
                .format(os.path.curdir))
    db.init(db_name)
    db.connect()
    yield db
    db.close()


def create_tables(drop=True):
    '''
    Creates the tables in the database. If `drop`=True, attempts to safely drop
    the tables before creating them.
    '''
    tbl_list = [Primer, Set, PrimerSet, _metadata]
    if drop:
        db.drop_tables(tbl_list, safe=True)
    db.create_tables(tbl_list, safe=True)


def check_create_tables(primer_db, skip_check=False):
    if os.path.isfile(primer_db) and not skip_check:
        warn(
            "This directory was already initialized as a workspace. "
            "Continuing will reset any stored primers and sets you may have found."
        )
        click.confirm("Are you sure you want to proceed?", abort=True)
    create_tables()


def check_version(swga_version=__version__):
    try:
        meta = _metadata.get()
        db_ver = semver.Version(meta.version)
    except pw.OperationalError:
        db_ver = "<NA>"
    ver = semver.Version(swga_version)
    spec = semver.Spec('=={}.{}'.format(ver.major, ver.minor))
    if db_ver not in spec:
        error(
            "This workspace was created with an incompatible version of swga.\n"
            "  Workspace version: {}\n"
            "  swga version:      {}\n"
            "Please re-initialize the workspace with `swga init` or use a "
            "different version of swga."
            .format(db_ver, ver),
            exception=False,
            wrap=False
        )


def add_primers(primers, chunksize=199, add_revcomp=True):
    if add_revcomp:
        def mkrevcomp(p):
            p2 = dict(**p)
            p2['seq'] = locate.revcomp(p['seq'])
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
        error("Invalid primers for set")
    if isinstance(primers, pw.SelectQuery):
        nprimers = primers.count()
    else:
        nprimers = len(primers)
    if nprimers == 0:
        error("Cannot have an empty set")
    _hash = hash(frozenset([p.seq for p in primers]))
    if Set.select(pw.fn.Count(Set._id)).where(Set._hash == _hash).scalar() > 0:
        return None
    s = Set.create(_id=_id, _hash=_hash, **kwargs)
    try:
        s.primers.add(primers)
    # If this is a Primers object, return the actual primers in it
    except AttributeError:
        s.primers.add(primers.primers)
    return s


def get_primers_for_set(set_id):
    set = Set.get(Set._id == set_id)
    return [p.seq for p in set.primers]


def get_primers_for_ids(pids):
    return list(Primer.select().where(Primer._id << pids).execute())


def get_metadata(db_name):
    m = _metadata.get()
    return meta(m.db_name, m.fg_file, m.bg_file, m.ex_file, m.fg_length, m.bg_length)

def set_metadata(db_name, version, fg_file, bg_file, ex_file, fg_length, bg_length):
    """Sets the metadata for a workspace. Clears any previous info."""
    _metadata.delete().execute()
    _metadata.insert(
        db_name=db_name,
        version=version,
        fg_file=fg_file,
        bg_file=bg_file,
        ex_file=ex_file,
        fg_length=fg_length,
        bg_length=bg_length
    ).execute()

