from swga.commands2.filter import Filter

def test_filter():
    cmd = Filter('filter', 'description')
    cmd.parse_args(['--input', 'infile',
                    '--output', 'outfile',
                    '--max_bg_binding', '5000',
                    '--num_primers', '200'])
    assert vars(cmd.args) == {'input': 'infile',
                        'output': 'outfile',
                        'max_bg_binding': 5000,
                        'num_primers': 200}
