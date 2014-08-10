PrimerSets
==========

Module to select sets of primers

Installation
-----------

This module depends on numpy > 1.6.1.

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

Default options for these commands are specified in `parameters.cfg`. All options can be overriden on the command line.

To begin, use [SelectiveWholeGenomeAmplification](https://github.com/mutantturkey/SelectiveWholeGenomeAmplification) to find a list of starting primers. This list should have one primer per line, with foreground binding numbers, background binding numbers, and weight ratio following, separated by spaces. For instance:
```
AATTGGCC 120000 1200 0.001
AATTCCGG 14000 2300 0.1642
```


## Theory

We need to select a compatible set of primers that binds to a target (foreground) genome with a higher frequency than a background genome for selective whole genome amplification. The primers need to be compatible with each other- i.e. not heterodimers- and have appropriate foreground and background binding frequencies. 

Given a list of primers, the goal of finding a group of _n_ compatible primers (where _n_ is the desired set size) is mappable to finding the largest [clique](https://en.wikipedia.org/wiki/Clique_(graph_theory)) in a graph.

A clique in graph theory is a set of nodes in an undirected graph that are all connected to one another. If we visualize all the primers as nodes on a graph, and connect each compatible pair of primers with an edge, we can see that a set of mutually compatible primers takes the form of a clique in this graph.

Unfortunately, this is a known NP-complete problem, so no known efficient solution exists. Exponential-time solutions do exist, however, and we will be using an implementation of the branch-and-bound algorithm called [cliquer](http://users.tkk.fi/~pat/cliquer.html) to find the desired primer sets.

Development will be focused on wrapping cliquer in Python to manage the workflow, devising a way of splitting the computational workload amongst multiple cores, and integrating with the rest of the SWGA pipeline.
