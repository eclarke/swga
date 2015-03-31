#### Messages used during `swga init`

welcome_message = """
## SWGA Initialization ---------------------------
This will set up a workspace for SWGA in the current directory ({CWD}).
"""

fg_message = """
Foreground genome: {fg_genome_fp}
  Length:  {fg_length} bp
  Records: {fg_nrecords}
"""

bg_message = """
Background genome: {bg_genome_fp}
  Length:  {bg_length} bp
"""

exclude_prompt = """
Do you want to add a FASTA file of sequences that will be used to exclude
primers? For instance, to avoid primers that bind to a mitochondrial genome,
you would add the path to that genome file. There can be multiple sequences in
the file, but only one file can be specified.
"""

finished_message = """
Done! 

A file called "{DEFAULT_PARAMETERS_FNAME}" has been placed in this directory. The
values given in this file will be used as defaults for the commands in SWGA. You
can modify these default values by editing this file in a plain-text editor such
as nano or TextEdit.
--------------------
"""
