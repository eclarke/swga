import click
from click._compat import filename_to_ui
import os, sys, stat
from pkg_resources import resource_string
from pyfaidx import Fasta
from swga.utils.options import load_swga_opts, cfg_from_opts

welcome_message = """
## SWGA Initialization ---------------------------

This will set up a workspace for SWGA in the current directory ({cwd}).

As part of this process, a default parameters file will be created that contains
some reasonable options for each part of the pipeline. This will be
named "{default_parameters_name}" and should be modified as appropriate for your
particular project. Please note that you should only with edit this file with a
plain-text editor like Notepad, TextEdit, or gedit; or with command-line tools
such as nano, vim, or emacs. 

The values specified in {default_parameters_name} will be used if the
corresponding values are not specified on the command line, or if the
"swga autopilot" command is used.
"""

fg_message = """
Foreground genome: {fg_genome}
  Length:  {fg_length} bp
  Records: {fg_nrecords}
"""
bg_message = """
Background genome: {bg_genome}
  Length:  {bg_length} bp
  Records: {bg_nrecords}
"""

finish_message = """Done!"""

@click.command()
@click.option("-f", "--fg_genome",
              type=click.Path(exists=True, resolve_path=True))
@click.option("-b", "--bg_genome",
              type=click.Path(exists=True, resolve_path=True))
def main(fg_genome, bg_genome):

    default_parameters_name = "default_parameters.cfg"
    cwd = os.getcwd()

    click.echo(click.style(welcome_message.format(**locals()),
                           fg = "blue"))
                                       
    if (not fg_genome):
        fg_genome = click.prompt("Enter path to foreground genome file, in " +
        "FASTA format", type=click.Path(exists=True, resolve_path=True))
    fg_length, fg_nrecords = fasta_stats(fg_genome)
    click.echo(click.style(fg_message.format(**locals()), fg="green"))

    if (not bg_genome):
        bg_genome = click.prompt("Enter path to background genome file, in " +
        "FASTA format", type=click.Path(exists=True, resolve_path=True))
    bg_length, bg_nrecords = fasta_stats(bg_genome)
    click.echo(click.style(bg_message.format(**locals()), fg="green"))

    opts = load_swga_opts()
    default_parameters = cfg_from_opts(opts)

    with open(os.path.join(cwd, default_parameters_name), "wb") as cfg_file:
        cfg_file.write(default_parameters.format(fg_genome_fp = fg_genome,
                                                 bg_genome_fp = bg_genome,
                                                 bg_length = bg_length))
    
    click.echo(finish_message)
              

def fasta_stats(fasta_fp):
    try:
        fasta = Fasta(fasta_fp)
        length = sum([len(r) for r in fasta])
        nrecords = len(fasta.keys())
        return length, nrecords
    except:
        click.echo(click.style("\nError reading %s: invalid FASTA format?" %
                               fasta_fp, fg = "red"))
        sys.exit(1)


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
