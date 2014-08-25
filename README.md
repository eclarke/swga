PrimerSets
==========

A pipeline for finding compatible sets of primers for selective whole genome amplification.

Requirements
-------------
- Linux or Mac OS X
- Python >= 2.7 (has not yet been tested on Python 3+)
- C compiler such as `gcc`


Installation
-----------

Download the most recent copy of the source [here](https://github.com/BrissonEEDS/PrimerSets/archive/master.zip). From a terminal:
```sh
unzip PrimerSets-master.zip
# Compiling the C extensions
cd PrimerSets-master/PrimerSets/lib
make
cd ..
# To install for all users:
sudo python setup.py install
# Alternatively, install for one user:
python setup.py install --user
```

To ensure everything worked correctly, you can run tests:
```sh
# from the downloaded PrimerSets-master/PrimerSets directory
python test/tests.py
```

After installation completes successfully, you should be able to call the various pipeline commands from any directory.

Configuration
-------------
All the options for the pipeline can be specified in a configuration file. An example config file is included in at `PrimerSets/parameters.cfg`- it is recommended to copy this file to your project's working directory and modify the values as needed.

By default, the script looks for this config file in the current working directory. However, if you want a global configuration file, you can set the `$swga_params` environmental variable and the scripts will find it. To keep this value persistent across sessions, add this to your .bashrc file or equivalent:
```sh
export swga_params=$HOME/path_to_cfg_file/parameters.cfg
```

Each option in the config file can be overridden on the command line.

# Quick start

Assuming you have a list of primers together with foreground and background binding rates and ratio in a file called `selected-mers`, a foreground genome at "fg-genome.fasta" and you have compiled the set_finder as above, the following commands will do the following:


```sh
# select at most 200 primers with fewer than 12000 bg binding sites
filter_primers -i selected-mers -o filtered_primers

# flatten the foreground genome for easier searching
fasta_flattener fg-genome.fasta > fg-genome.fasta.flattened

# find locations of filtered primers in foreground genome
find_fg_locations -v -i filtered_primers --fg_genome fg-genome.fasta.flattened

# remove heterodimers and output compatibility graph
mk_primer_graph -i filtered_primers > primer_graph

# use compatibility graph to find sets btwn 2 and 7 primers,
# then filter for sets with less than 36000 bp btwn foreground binding sites
find_sets -i primer_graph | process_sets -o valid_sets.txt
```
Right now the `valid_sets` file is in a space-delimited format where the first col is set standard deviation, second is the max foreground genome binding distance, and remaining numbers are the indexes of the primers in `filtered_primers`.

If you wanted to be really cool, you can do all of this in only two lines:
```sh
fasta_flattener fg-genome.fasta > fg-genome.fasta.flattened
filter_primers -i selected-mers | find_fg_locations -p | mk_primer_graph | find_sets | \
    process_sets -o valid_sets.txt
```
since most of these commands will default to stdin/stdout if input and output are not specified.

# Usage
`PrimerSets` has 5 pipeline scripts. In order of usage:
- `filter_primers`: filters primers from input that match criteria
- `find_fg_locations`: finds locations of primers on the foreground genome and stores results
- `mk_primer_graph`: creates the primer compatibility graph before finding sets (see [theory](#Theory))
- `find_sets`: finds initial sets of compatible primers using branch-and-bound algorithm
- `process_sets`: filters sets that match criteria

Default options for these commands are specified in `parameters.cfg`. All options can be overridden on the command line; type `<script> -h` for options. Since most scripts accept input from stdin if not otherwise specified, running them without arguments will not show a help message.

To begin, use [SelectiveWholeGenomeAmplification](https://github.com/mutantturkey/SelectiveWholeGenomeAmplification) to find a list of starting primers. This list should have one primer per line, with foreground binding numbers, background binding numbers, and weight ratio following, separated by spaces or tabs. For instance:
```
AATTGGCC 120000 1200 0.001
AATTCCGG 14000 2300 0.1642
```

1.  The first step is to run `filter_primers` on this list. This removes any primers that bind to the background
    genome more than the number of times specified either on the command line or in the config file. It then
    reorders them according to the fourth column, which is usually a modified ratio of the number of
    foreground/background binding sites. Direct the output to a file with the `--output` flag.
    ```sh
    filter_primers -i selected-mers -o filtered_primers
    ```

2. The next step involves finding the binding sites of these filtered primers on the foreground genome. For the
    sake of easier coding, this module searches through a "flat" FASTA file, which is basically a FASTA file
    with all the header text and newlines removed. Use the `fasta_flattener` tool to convert your foreground
    genome, and update `parameters.cfg` with the flattened FASTA file's location.
    ```sh
    fasta_flattener fg-genome.fasta > fg-genome.fasta.flattened
    ```

3. To actually find the locations of the primers, run `find_fg_locations -i filtered_primers` where
    `filtered_primers` is the output file from step 1. For particularly large foreground genomes, you can add
    `-v` to the command to track progress. To chain this command with other scripts, add the `-p` or
    `--passthrough` flag, which echos the input to stdout (see the two-line example in the
    [Quick Start](#Quick start)).
    ```sh
    find_fg_locations -v -i filtered_primers --fg_genome fg-genome.fasta.flattened
    ```

4. We can now create the primer compatibility graph. This graph makes each primer a vertex in the graph and adds
    an edge between two primer vertices if they are not heterodimers (according to the thresholds specified).
    This is done with the `mk_primer_graph` command, and is usually extremely fast. As with all the other
    scripts, you can save the output with the `-o` or `--output` flag, or simply pipe the output directly to
    `find_sets`.
    ```sh
    find_fg_locations -v -i filtered_primers --fg_genome fg-genome.fasta.flattened
    ```

5. We now analyze the primer graph to find groups of mutually compatible primers using the `find_sets` script.  
    It takes the primer compatibility graph and attempts to find completely connected subgraphs (cliques) where
    the cumulative background genome binding frequency of the primers in the clique is lower than a certain
    amount. Since no efficient algorithm exists for this step, we use a highly optimized branch-and-bound
    algorithm implemented in C called [cliquer](http://users.tkk.fi/~pat/cliquer.html). `find_sets` is a wrapper
    around a modified version of `cliquer`. For more on this step, see the [Theory](#theory) section below.
    ```sh
    find_fg_locations -v -i filtered_primers --fg_genome fg-genome.fasta.flattened
    ```

6. The `process_sets` command takes the output from `find_sets` and runs further statistics on the sets, only
    outputting a set if the primers in the set have a max distance on the foreground genome below a certain
    threshold. It also currently calculates the standard deviation of the binding distances on the foreground
    genome. The output of this command has the format per row of
    `[stdev] [max_dist] [primer seq 1] [primer seq 2] [...]`.
    ```sh
    find_fg_locations -v -i filtered_primers --fg_genome fg-genome.fasta.flattened
    ```

There are often billions of valid sets possible, and it would be extremely time-consuming to find all of them. Therefore, `process_sets` can read the output from `find_sets` live and will exit after it finds `max_sets` valid primer sets. Alternatively, one could let the `find_sets` write to a file and find all sets, then pass it to `process_sets` and test different parameters.


## Theory

We need to select a compatible set of primers that binds to a target (foreground) genome with a higher frequency than a background genome for selective whole genome amplification. The primers need to be compatible with each other- i.e. not heterodimers- and have appropriate foreground and background binding frequencies.

Given a list of primers, the goal of finding a group of _n_ compatible primers (where _n_ is the desired set size) is mappable to finding the largest [clique](https://en.wikipedia.org/wiki/Clique_(graph_theory)) in a graph.

A clique in graph theory is a set of nodes in an undirected graph that are all connected to one another. If we visualize all the primers as nodes on a graph, and connect each compatible pair of primers with an edge, we can see that a set of mutually compatible primers takes the form of a clique in this graph.

Unfortunately, this is a known NP-complete problem, so no known efficient solution exists. Exponential-time solutions do exist, however, and we will be using an implementation of the branch-and-bound algorithm called [cliquer](http://users.tkk.fi/~pat/cliquer.html) to find the desired primer sets.

# Contributors and Acknowledgements

## Coding:
  - [Erik Clarke](https://github.com/eclarke) -- main implementation
  - [Calvin Morrison](https://github.com/mutantturkey) -- heterodimer comparison
