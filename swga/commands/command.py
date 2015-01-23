import swga

class Command:


    def __init__(self, name, description, cfg_file):
        self.parser = swga.basic_cmd_parser(description, name, cfg_file)


    def main(self, argv):
        self.args = self.parser.parse_args(argv)
    

    def run(self, args):
        pass
