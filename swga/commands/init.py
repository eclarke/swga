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
from click._compat import filename_to_ui
import os
import stat
import swga
from pyfaidx import Fasta
from swga.utils import (resources, options)
#from swga.utils.resources import get_swga_opts
#from swga.utils.options import cfg_from_opts
from swga.data.messages import (
    welcome_message,
    fg_message,
    bg_message,
    exclude_prompt,
    finished_message)

SWGA_VERSION = swga.__version__

DEFAULT_PARAMETERS_FNAME = "parameters.cfg"

# Calculate default binding frequencies for foreground and background
# Fg is based on a binding rate of 1/100000 bp/binding site
# Bg is based on a binding rate of 1/150000 bp/binding site
MIN_FG_RATE = 0.00001
MAX_BG_RATE = 0.0000067


@click.command()
@click.option(
    "-f", "--fg_genome_fp",
    type=click.Path(exists=True, resolve_path=True))
@click.option(
    "-b", "--bg_genome_fp",
    type=click.Path(exists=True, resolve_path=True))
@click.option(
    "-e", "--exclude_fp",
    type=click.Path(exists=True, resolve_path=True))
def main(fg_genome_fp, bg_genome_fp, exclude_fp):
    CWD = os.getcwd()

    # 01. Display welcome message
    click.secho(
        welcome_message.format(
            CWD=CWD, SWGA_VERSION=SWGA_VERSION), fg="blue")

    # 02. Prompt for the foreground genome, if not already specified
    if (not fg_genome_fp):
        fg_genome_fp = click.prompt(
            "Enter path to foreground genome file, in FASTA format",
            type=click.Path(exists=True, resolve_path=True))
    fg_length, fg_nrecords = fasta_stats(fg_genome_fp)
    click.secho(fg_message.format(**locals()), fg="green")

    # 03. Prompt for background genome, if not already specified
    if (not bg_genome_fp):
        bg_genome_fp = click.prompt(
            "Enter path to background genome file, in FASTA format",
            type=click.Path(exists=True, resolve_path=True))
    bg_length = fasta_len_quick(bg_genome_fp)
    click.secho(bg_message.format(**locals()), fg="green")

    # 04. Prompt for a file containing sequences to exclude, if not already
    # given
    if (not exclude_fp):
        if click.confirm(exclude_prompt):
            exclude_fp = click.prompt(
                "Enter path to exclusionary sequence(s), in FASTA format",
                type=click.Path(exists=True, resolve_path=True))

    if exclude_fp:
        exclude_fp_message = click.style(
            "Exclusionary sequences file: {}".format(exclude_fp), fg="red")
    else:
        exclude_fp = ""
        exclude_fp_message = click.style(
            "No exclusionary sequences file specified.", fg="green")

    click.echo(exclude_fp_message)

    # 05. Build and populate the parameters file
    opts = resources.get_swga_opts()
    default_parameters = options.cfg_from_opts(opts)
    cfg_fp = os.path.join(CWD, DEFAULT_PARAMETERS_FNAME)
    min_fg_bind = int(MIN_FG_RATE * float(fg_length))
    max_bg_bind = int(MAX_BG_RATE * float(bg_length))

    # 06. Write parameters file
    if os.path.isfile(cfg_fp):
        click.confirm(
            "Existing file `%s` will be overwritten. Continue?"
            % DEFAULT_PARAMETERS_FNAME, abort=True)

    with open(cfg_fp, "wb") as cfg_file:
        cfg_file.write(default_parameters.format(SWGA_VERSION=SWGA_VERSION, **locals()))

    # Done!
    click.secho(finished_message.format(
        DEFAULT_PARAMETERS_FNAME=DEFAULT_PARAMETERS_FNAME), fg="green")


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
    creates a FASTA index (.fai) for later searching. May be slow for very large
    files.
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


if __name__ == "__main__":
    main()
