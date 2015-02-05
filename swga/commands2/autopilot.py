
import swga
import pipes
from swga.commands2 import Command


def main(argv, cfg_file):
    cmd = Command('autopilot', cfg_file=cfg_file)
    cmd.parse_args(argv)
    autopilot(**cmd.args)


def autopilot():
    ap = pipes.Template()
    ap.prepend('swga find_sets', '--')
    ap.append('swga score 1> $OUT', '-f')
    f = "set_scores.txt"
    ap.open(f, 'w')
    
