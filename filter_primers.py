from PrimerSets import *
from signal import signal, SIGPIPE, SIG_DFL
import subprocess
from os.path import isfile


def main():
    opts = read_config_file(default_config_file)
    if not isfile(default_config_file):
        opts.add_section('primer_filters')
        opts.set('primer_filters', 'max_bg_binding', '0')
        opts.set('primer_filters', 'num_primers', '0')
        sys.stderr.write(opts_errstr)
    
    argparser = ArgumentParser(description = '''Reads in a list of
    primers, removing those that bind to the background genome over a
    threshold, and ordering the remaining by the fg/bg binding
    ratio. Defaults values for threshold and number of primers in the
    config file.''') 

    argparser.add_argument('primer_file', action='store', help='''List of
    primers (first col), number of foreground binding sites (second
    col), number of background binding sites (third col), and fg/bg
    ratio (fourth col)''')
    
    max_bg_binding_default = opts.getint('primer_filters', 'max_bg_binding')
    argparser.add_argument('-m', '--max_bg_binding', help="""Max times
    a primer is allowed to bind to the background genome.""", 
                           default=max_bg_binding_default)

    num_primers_default = opts.getint('primer_filters', 'num_primers')
    argparser.add_argument('-n', '--num_primers', help="""Max number
    of primers to use (may be less).""",
                           default=num_primers_default) 

    args = vars(argparser.parse_args())

    primer_file = args['primer_file']
    max_bg_binding = args['max_bg_binding']
    num_primers = args['num_primers']
    
    # Sorts the input by the third col (bg binding), trims if value
    # more than max_bg_binding, sorts results by fourth col (ratio),
    # and prints the first $num_primers of the result.
    # This could be done in python but the succinctness is alluring.
    command = """sort -t ' ' -n -k 3 < {} | awk '{{if ($3 < {}) print
    $0}}' | sort -t ' ' -n -r -k 4 | head -n {}""".format(primer_file,
    max_bg_binding, num_primers) 

    # the preexec_fn is in there to prevent sort complaining about a
    # broken pipe
    subprocess.call(command, shell=True,
                    preexec_fn = lambda: signal(SIGPIPE, SIG_DFL))
    

if __name__ == "__main__":
    main()
