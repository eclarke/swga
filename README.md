PrimerSets
==========

Find compatible sets of primers for selective whole genome amplification.

Installation
-----------

This program has no external dependencies besides a C compiler such as `gcc`.

Download the most recent copy of the source [here](https://github.com/BrissonEEDS/PrimerSets/archive/master.zip). From a terminal:
```sh
unzip PrimerSets-master.zip
cd PrimerSets-master/swga/cliquer
make all
make install
```
# Quick start

Assuming you have a list of primers together with foreground and background binding rates and ratio in a file called `selected-mers`, a foreground genome at "fg-genome.fasta" and you have compiled the set_finder as above, the following commands will do the following:


```sh
cd swga

# select at most 200 primers with fewer than 12000 bg binding sites
./swga.py filter_primers selected-mers > filtered_primers

# flatten the foreground genome for easier searching
utils/fasta_flattener.sh fg-genome.fasta > fg-genome.fasta.flattened

# find locations of filtered primers in foreground genome
./swga.py fg_locations -v -i filtered_primers \
  --fg_genome fg-genome.fasta.flattened

# remove heterodimers and output compatibility graph
./swga.py make_graph -i filtered_primers > primer_graph

# use compatibility graph to find sets btwn 2 and 7 primers, 
# then filter for sets with less than 36000 bp btwn foreground binding sites 
./swga.py find_sets -i primer_graph | ./swga.py process_sets -o valid_sets.txt
```
Right now the `valid_sets` file is in a space-delimited format where the first col is set standard deviation, second is the max foreground genome binding distance, and remaining numbers are the indexes of the primers in `filtered_primers`.

If you wanted to be really cool, you can do all of this in only two lines:
```shell
utils/fasta_flattener.sh fg-genome.fasta > fg-genome.fasta.flattened
./swga.py filter_primers selected-mers | ./swga.py fg_locations -v \
    | ./swga.py make_graph | ./swga.py find_sets | ./swga.py process_sets -o valid_sets.txt
```
since most of these commands will default to stdin/stdout if input and output are not specified.

# Usage
From the `swga` directory, run `./swga.py command`.

`swga.py` has 5 command modes:
- `filter_primers`: filters primers from input that match criteria
- `fg_locations`: finds locations of primers on the foreground genome and stores results
- `make_graph`: creates the primer compatibility graph before finding sets (see [theory](#Theory))
- `find_sets`: finds initial sets of compatible primers using branch-and-bound algorithm
- `process_sets`: filters sets that match criteria

Default options for these commands are specified in `parameters.cfg`. All options can be overriden on the command line; type `swga.py <command> -h` for options.

To begin, use [SelectiveWholeGenomeAmplification](https://github.com/mutantturkey/SelectiveWholeGenomeAmplification) to find a list of starting primers. This list should have one primer per line, with foreground binding numbers, background binding numbers, and weight ratio following, separated by spaces. For instance:
```
AATTGGCC 120000 1200 0.001
AATTCCGG 14000 2300 0.1642
```

The first step is to run `swga.py filter_primers` on this list. This removes any primers that bind to the background genome more than the number of times specified either on the command line or in the config file. It then reorders them according to the fourth column, which is usually a modified ratio of the number of foreground/background binding sites. The command usually outputs straight to standard out, so to save its output, simply redirect it with `> output_file`.

The next step involves finding the binding sites of these filtered primers on the foreground genome. For the sake of easier coding, this module searches through a "flat" FASTA file, which is basically a FASTA file with all the header text and newlines removed. Use the `utils/fasta_flattener.sh` tool to convert your foreground genome, and update `parameters.cfg` with the flattened FASTA file's location. Next, run `swga.py fg_locations -i filtered_primers` where `filtered_primers` is the file from step 1. Normally this command echos its input to stdout so you can chain it together with other commands. To avoid this behavior, simply add `--no_passthrough` to the command. For particularly large foreground genomes, you can add `-v` to the command to track progress.

At this point, we can now create the primer compatibility graph. This graph contains each primer as a vertex and adds an edge between two primer vertices if they are not heterodimers (according to the thresholds specified). This is done with the `swga.py make_graph` command, and is usually extremely fast, so there is no need to save its output. Instead, it is equally easy to chain it together with the `swga.py find_sets` and `swga.py process_sets` commands, like so:
```sh
swga.py make_graph -i filtered_primers | swga.py find_sets | swga.py process_sets -o valid_sets.txt
```

The `find_sets` command is the heart of this module. It takes the primer compatibility graph and attempts to find completely connected subgraphs (cliques) where the cumulative background genome binding frequency of the primers in the clique is lower than a certain amount. See [theory](#theory) below for more explanation.

The `process_sets` command takes the output from `find_sets` and runs further statistics on the sets, only outputting a set if the primers in the set have a max distance on the foreground genome below a certain threshold. It also currently calculates the standard deviation of the binding distances on the foreground genome. The output of this command has the format per row of `[stdev] [max_dist] [primer seq 1] [primer seq 2] [...]`.

There are often billions of valid sets possible, and it would be extremely time-consuming to find all of them. Therefore, `process_sets` can read the output from `find_sets` live and will exit after it finds `max_sets` valid primer sets. Alternatively, one could let the `find_sets` write to a file and find all sets, then pass it to `process_sets` and test different parameters.


## Theory

We need to select a compatible set of primers that binds to a target (foreground) genome with a higher frequency than a background genome for selective whole genome amplification. The primers need to be compatible with each other- i.e. not heterodimers- and have appropriate foreground and background binding frequencies. 

Given a list of primers, the goal of finding a group of _n_ compatible primers (where _n_ is the desired set size) is mappable to finding the largest [clique](https://en.wikipedia.org/wiki/Clique_(graph_theory)) in a graph.

A clique in graph theory is a set of nodes in an undirected graph that are all connected to one another. If we visualize all the primers as nodes on a graph, and connect each compatible pair of primers with an edge, we can see that a set of mutually compatible primers takes the form of a clique in this graph. 

Unfortunately, this is a known NP-complete problem, so no known efficient solution exists. Exponential-time solutions do exist, however, and we will be using an implementation of the branch-and-bound algorithm called [cliquer](http://users.tkk.fi/~pat/cliquer.html) to find the desired primer sets.

# Contributors

### Coding:
  - [Erik Clarke](https://github.com/eclarke) -- main implementation
  - [Calvin Morrison](https://github.com/mutantturkey) -- heterodimer comparison
