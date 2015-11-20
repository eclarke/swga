# -*- coding: utf-8 -*-
"""
`swga init` is a special command that initializes a directory before being used
with SWGA. In short, this reads the foreground and background FASTA files to
get some basic stats (number of records, overall number of bases) and
populates the generated parameters.cfg file with reasonable defaults based on
these stats.

For this reason, this command does not use the framework that other commands
use (since it creates the parameters file).
"""

import click
import subprocess
import argparse
from click._compat import filename_to_ui
import os
import stat
from swga import (
    DEFAULT_CFG_FNAME,
    DEFAULT_DB_FNAME,
    __version__
)
from swga.commands import create_config_file
import swga.database as database
from pyfaidx import Fasta

version = __version__

# Default binding frequencies for foreground and background
# Fg is based on a binding rate of 1/100000 bp/binding site
# Bg is based on a binding rate of 1/150000 bp/binding site
MIN_FG_RATE = 0.00001
MAX_BG_RATE = 0.0000067


setup_msg = '''
swga v{version} - interactive setup
-------------------------------
This will set up a new swga workspace in the current directory ({CWD}).
'''

fg_prompt = '''Enter path to foreground FASTA file'''
fg_parsing = 'Checking {fg_genome_fp}...'
fg_msg = '''\
Foreground filepath: {fg_genome_fp}
  Length:  {fg_length} bp
  Records: {fg_nrecords}
'''

bg_prompt = '''Enter path to background FASTA file'''
bg_parsing = 'Checking {bg_genome_fp}...'
bg_msg = '''\
Background filepath: {bg_genome_fp}
  Length:  {bg_length} bp
'''

excl_prompt = '''\
Do you want to add a FASTA file of sequences that will be used to exclude
primers? For instance, to avoid primers that bind to a mitochondrial genome,
you would add the path to that genome file. There can be multiple sequences in
the file, but only one file can be specified.
'''
excl_prompt2 = '''Enter path to exclusionary sequence(s), in FASTA format'''
excl_msg = '''Exclusionary sequences file: {}'''
excl_msg_no_file = '''No exclusionary sequences file specified.'''

overwrite_cfg_prompt = '''Existing file `{}' will be overwritten. Continue?'''

fin_msg = '''\
Done!

Created pre-filled config file `{}'. Workspace data: `{}'

This file has been pre-filled with reasonable defaults. However, you will
probably want to modify them to suit your needs. You can override the values
in this file both by passing arguments to each command or by modifying this
file in a plain-text editor such as TextEdit or `nano'.
--------------------
'''


def main(argv=None):

    parser = argparse.ArgumentParser(
        description="Initialize the current directory as an swga workspace.")
    parser.add_argument(
        '-f', '--fg_genome_fp',
        help='Path to foreground genome/sequences (in FASTA format)',
        type=argparse.FileType('r'))
    parser.add_argument(
        '-b', '--bg_genome_fp',
        help='Path to background genome/sequences (in FASTA format)',
        type=argparse.FileType('r'))
    parser.add_argument(
        '-e', '--exclude_fp',
        help='Path to sequences to exclude from analysis (in FASTA format)',
        type=argparse.FileType('r'))
    parser.add_argument(
        '--force',
        help='Overwrite any existing data in the directory without prompts',
        action='store_true')

    args = parser.parse_args(argv) if argv else parser.parse_args()

    # Convert to local variables
    fg_genome_fp = args.fg_genome_fp
    bg_genome_fp = args.bg_genome_fp
    exclude_fp = args.exclude_fp

    # Sanity check our default files for read/write
    cwd = os.path.abspath(os.curdir)
    assert os.access(cwd, os.W_OK | os.X_OK)
    db_fp = os.path.abspath(os.path.join(cwd, DEFAULT_DB_FNAME))
    # assert os.access(db_fp, os.W_OK | os.X_OK)
    cfg_fp = os.path.abspath(DEFAULT_CFG_FNAME)
    # assert os.access(cfg_fp, os.W_OK | os.X_OK)

    # 01. Display welcome message
    click.secho(setup_msg.format(version=version, CWD=cwd), fg="blue")

    # 02. Prompt for the foreground genome, if not already specified
    if (not fg_genome_fp):
        fg_genome_fp = click.prompt(
            fg_prompt, type=click.Path(exists=True, resolve_path=True))
    else:
        fg_genome_fp = os.path.abspath(fg_genome_fp.name)
    click.secho(fg_parsing.format(**locals()), fg='blue')
    fg_length, fg_nrecords = fasta_stats(fg_genome_fp)
    click.secho(fg_msg.format(**locals()), fg="green")

    # 03. Prompt for background genome, if not already specified
    if (not bg_genome_fp):
        bg_genome_fp = click.prompt(
            bg_prompt, type=click.Path(exists=True, resolve_path=True))
    else:
        bg_genome_fp = os.path.abspath(bg_genome_fp.name)
    click.secho(bg_parsing.format(**locals()), fg='blue')
    bg_length = fasta_len_quick(bg_genome_fp)
    click.secho(bg_msg.format(**locals()), fg="green")

    # 04. Prompt for a file containing sequences to exclude, if not already
    #     given
    if (not exclude_fp):
        if click.confirm(excl_prompt):
            exclude_fp = click.prompt(
                excl_prompt2, type=click.Path(exists=True, resolve_path=True))
    else:
        exclude_fp = os.path.abspath(exclude_fp.name)

    if exclude_fp:
        click.secho(excl_msg.format(exclude_fp), fg="green")
    else:
        exclude_fp = ""
        click.secho(excl_msg_no_file, fg="green")

    # 05. Build and populate the config file
    default_parameters = create_config_file()
    #cfg_fp = os.path.join(CWD, DEFAULT_CFG_FNAME)
    min_fg_bind = int(MIN_FG_RATE * float(fg_length))
    max_bg_bind = int(MAX_BG_RATE * float(bg_length))

    # 06. Write config file
    if not args.force and os.path.isfile(cfg_fp):
        click.confirm(overwrite_cfg_prompt.format(DEFAULT_CFG_FNAME),
                      abort=True)
    with open(cfg_fp, "wb") as cfg_file:
        cfg_file.write(default_parameters.format(**locals()))

    # 07. Initialize the database
    database.init_db(db_fp, create_if_missing=True)
    database.check_create_tables(db_fp, skip_check=args.force)
    database.set_metadata(
        db_fp,
        version,
        fg_genome_fp,
        bg_genome_fp,
        exclude_fp,
        fg_length,
        bg_length)

    # Done!
    click.secho(fin_msg.format(cfg_fp, db_fp), fg="green")


def fasta_len_quick(fasta_fp):
    """
    Fast way to get the number of bases in a FASTA file, excluding headers.
    """
    try:
        length = subprocess.check_output(
            "grep -v '>' {} | wc -c".format(fasta_fp), shell=True)
        length = int(length.strip())
        return length
    except subprocess.CalledProcessError:
        click.echo(click.style(
            "\nError getting length of FASTA file {}".format(fasta_fp),
            fg="red"))
        raise


def fasta_stats(fasta_fp):
    """
    Retrieves the number of bases and number of records in a FASTA file. Also
    creates a FASTA index (.fai) for later searching. May be slow for very
    large files.
    """
    # pyfaidx can't handle blank lines within records, so we have to check :(
    check_empty_lines(fasta_fp)
    try:
        fasta = Fasta(fasta_fp)
        length = fasta_len_quick(fasta_fp)
        nrecords = len(fasta.keys())
        return length, nrecords
    except:
        click.secho(
            "\nError reading %s: invalid FASTA format?" % fasta_fp, fg="red")
        raise


def check_empty_lines(fasta_fp):
    '''
    The pyfaidx module has issues with blank lines in Fasta files. This checks
    for their presence and offers to remove them in-place.
    '''
    basename = os.path.basename(fasta_fp)
    try:
        check = subprocess.check_output(
            "grep -n '^$' {}".format(fasta_fp), shell=True)
    except subprocess.CalledProcessError:
        check = None
    if check:
        click.confirm(
            "`{}` has blank lines, which can mess up our FASTA indexer."
            " Is it okay to remove these lines from the file?"
            .format(basename), abort=True)
        # sed -ie gets around the fact that BSD sed -i requires a backup extension
        # with sed, but we don't want to create a backup. Ensures linux/osx
        # compatibility.
        subprocess.check_call(
            "sed -ie '/^$/d' \"{}\"".format(fasta_fp), shell=True)

##
# The monkeypatching going on here forces the click framework to strip
# whitespace from the ends of filename strings passed to it. This allows users
# to drag-and-drop files on the Mac OS X terminal (which appends a trailing
# space); without this, click says that the file doesn't exist.
##


def monkeypatch_method(cls):
    def decorator(func):
        setattr(cls, func.__name__, func)
        return func
    return decorator


@monkeypatch_method(click.Path)
def convert(self, value, param, ctx):
    value = value.strip()
    rv = value
    if self.resolve_path:
        rv = os.path.realpath(rv)

    try:
        st = os.stat(rv)
    except OSError:
        if not self.exists:
            return rv
        self.fail('%s "%s" does not exist.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)

    if not self.file_okay and stat.S_ISREG(st.st_mode):
        self.fail('%s "%s" is a file.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if not self.dir_okay and stat.S_ISDIR(st.st_mode):
        self.fail('%s "%s" is a directory.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if self.writable and not os.access(value, os.W_OK):
        self.fail('%s "%s" is not writable.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if self.readable and not os.access(value, os.R_OK):
        self.fail('%s "%s" is not readable.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)

    return rv
