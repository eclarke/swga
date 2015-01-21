from swga.clint.textui import puts, colored, indent
from swga.clint.arguments import Args

args = Args()

message = """
Install succeeded!
{user_message}
To get started, run:
{swga_cmd}

For help, see the docs at https://github.com/eclarke/swga/wiki
"""

def main():
    cmd_prefix = args.get(0) if args.get(0) else ""
    swga_cmd = colored.green(cmd_prefix+"swga")
    user_message = ""
    if cmd_prefix:
        user_message = colored.blue("\nWe recommend putting\n {prefix}\nin your path.\n").format(prefix=colored.magenta(cmd_prefix))
    puts(colored.blue('#'*70))
    with indent(2, quote=colored.blue("#")):
        puts(colored.blue(message).format(swga_cmd = swga_cmd,
                                          user_message = user_message))
    puts(colored.blue('#'*70))

if __name__ == "__main__":
    main()
