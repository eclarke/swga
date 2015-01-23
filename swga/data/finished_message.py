from swga.clint.textui import puts, colored, indent
from swga.clint.arguments import Args

args = Args()

message = """
Path to SWGA:
{swga_cmd}

Either use the full path each time you want to run SWGA, or add 
{prefix} to your path by adding
\texport PATH=$PATH:{prefix}
to your ~/.bashrc file.

For help, see the docs at https://github.com/eclarke/swga/wiki
"""

def main():
    cmd_prefix = args.get(0) if args.get(0) else ""
    swga_cmd = colored.green(cmd_prefix+"swga")
    puts(colored.blue('#'*70))
    with indent(2, quote=colored.blue("#")):
        puts(message.format(swga_cmd = swga_cmd,
                            prefix = colored.magenta(cmd_prefix)))
    puts(colored.blue('#'*70))

if __name__ == "__main__":
    main()
