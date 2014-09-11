PrimerSets
==========

A pipeline for finding compatible sets of primers for selective whole genome amplification.

Requirements
-------------
- Linux or Mac OS X
- Python >= 2.7 (has not yet been tested on Python 3+)
- pip >= 1.5.6
- setuptools >= 5.7
- C compiler such as `gcc`


Installation
-----------

The easiest way to download and keep the pipeline up-to-date is to install
an editable copy of it from the repository. Just clone a copy and run `make`:
```sh
git clone https://github.com/BrissonEEDS/PrimerSets
cd PrimerSets
make
```
By default, the installation uses `~/.swga` as the location of `SWGAHOME`, which is where it keeps needed resources. To use a different directory for `SWGAHOME`, alter the above instructions so that the make command looks like this:
```sh
make SWGAHOME=/my/other/directory
```

After installation, if everything worked, you should see a message that looks something like this (from my Mac- your message may look different):

```
------------------------------------
Install succeeded!
Before using SWGA, run the following commands or put them in your
shell config file (such as .bashrc, .profile, or .bash_profile).
If not placed in your shell config file, these commands must be run each
time before using SWGA:

export SWGAHOME=/Users/erik/.swga
export PATH=$PATH:/Users/erik/Library/Python/2.7/bin
------------------------------------
```

The two statements starting with `export` should be added to your shell config file. On Linux systems, this may be ~/.profile or ~/.bashrc. On a Mac, this may be ~/.bash_profile. To add these statements, open the correct file for your system using an editor such as `nano` and paste the two statements at the bottom of that file, then save and quit.

Now, open a new terminal window and verify that everything worked:
```sh
# This should have the new directory at the end
echo $PATH
# This should return the directory given in the second message for SWGAHOME
echo $SWGAHOME
# Try running swga!
swga
```

If everything worked, congrats! Time to set up your configuration file.

Configuration
-------------

All the options for the pipeline can be specified in a configuration file. The pipeline first looks for a config file specified by the `--config` flag on the command line. If unspecified, it looks for a file called `parameters.cfg` in the local directory.

*Note*: Because there are many options for the pipeline, it is best practice (and now required) to specify a config file (either the default in the working directory, or by using `--config`).

Each option in the config file can be overridden on the command line.

# Quick start

Assuming you have a list of primers together with foreground and background binding rates and ratio in a file called `selected-mers`, a foreground genome at "fg-genome.fasta" and you have compiled the set_finder as above, the following commands will do the following:


```sh
# select at most 200 primers with fewer than 12000 bg binding sites
swga filter --input selected-mers --output filtered_primers

# flatten the foreground genome for easier searching
swga flatten --input fg-genome.fasta --output fg-genome.fasta.flattened

# find locations of filtered primers in foreground genome
swga locate --input filtered_primers --genome fg-genome.fasta.flattened

# remove heterodimers and output compatibility graph
swga mkgraph --input filtered_primers --output primer_graph.gr

# use compatibility graph to find sets btwn 2 and 7 primers,
# then filter for sets with less than 36000 bp btwn foreground binding sites
swga sets --input primer_graph.gr | swga score --output valid_sets.txt
```
Right now the `valid_sets` file is in a tab- or space-delimited format where the first col is set standard deviation, second is the max foreground genome binding distance, and remaining numbers are the indexes of the primers in `filtered_primers`.

If you wanted to be really cool, you can do all of this in only two lines:
```sh
swga flatten --input fg-genome.fasta --output fg-genome.fasta.flattened
swga filter -i selected-mers | swga locate -p | swga mkgraph | swga sets | \
    swga score -o valid_sets.txt
```
since most of these commands will default to stdin/stdout if input and output are not specified.

# Usage
The commands for SWGA are as follows:
- `count`: not yet implemented
- `flatten`: reduce a fasta file to one line (faster searching)
- `filter`: filters primers from input that match criteria
- `locate`: finds locations of primers on the input genome and stores results
- `mkgraph`: creates the primer compatibility graph before finding sets (see [theory](#Theory))
- `sets`: finds initial sets of compatible primers using branch-and-bound algorithm
- `score`: filters and scores sets that match criteria

Default options for these commands are specified in `parameters.cfg`. All options can be overridden on the command line; type `swga <command> -h` for options. Since most scripts accept input from stdin if not otherwise specified, running them without arguments will not show a help message.

To begin, use [SelectiveWholeGenomeAmplification](https://github.com/mutantturkey/SelectiveWholeGenomeAmplification) to find a list of starting primers. This list should have one primer per line, with foreground binding numbers, background binding numbers, and weight ratio following, separated by spaces or tabs. For instance:
```
AATTGGCC 120000 1200 0.001
AATTCCGG 14000 2300 0.1642
```

1.  The first step is to run `swga filter` on this list. This removes any primers that bind to the background
    genome more than the number of times specified either on the command line or in the config file. It then
    reorders them according to the fourth column, which is usually a modified ratio of the number of
    foreground/background binding sites. Direct the output to a file with the `--output` flag.
    ```sh
    swga filter -i selected-mers -o filtered_primers
    ```

2. The next step involves finding the binding sites of these filtered primers on the foreground genome. For the
    sake of easier coding, this module searches through a "flat" FASTA file, which is basically a FASTA file
    with all the header text and newlines removed. Use the `swga flatten` command to convert your foreground
    genome, and update `parameters.cfg` with the flattened FASTA file's location.
    ```sh
    swga flatten -i fg-genome.fasta -o fg-genome.fasta.flattened
    ```

3. To actually find the locations of the primers, run `swga locate -i filtered_primers` where
    `filtered_primers` is the output file from step 1. For particularly large foreground genomes, you can add
    `-v` to the command to track progress. To chain this command with other scripts, add the `-p` or
    `--passthrough` flag, which echos the input to stdout (see the two-line example in the
    [Quick Start](#Quick start)).
    ```sh
    swga locate -i filtered_primers --genome fg-genome.fasta.flattened
    ```

4. We can now create the primer compatibility graph. This graph makes each primer a vertex in the graph and adds
    an edge between two primer vertices if they are not heterodimers (according to the thresholds specified).
    This is done with the `swga mkgraph` command, and is usually extremely fast. As with all the other
    scripts, you can save the output with the `-o` or `--output` flag, or simply pipe the output directly to
    `swga sets`.
    ```sh
    swga mkgraph -i filtered_primers -o primer_graph.gr
    ```

5. We now analyze the primer graph to find groups of mutually compatible primers using the `swga sets` script.  
    It takes the primer compatibility graph and attempts to find completely connected subgraphs (cliques) where
    the cumulative background genome binding frequency of the primers in the clique is lower than a certain
    amount. Since no efficient algorithm exists for this step, we use a highly optimized branch-and-bound
    algorithm implemented in C called [cliquer](http://users.tkk.fi/~pat/cliquer.html). `find_sets` is a wrapper
    around a modified version of `cliquer`. For more on this step, see the [Theory](#theory) section below.
    ```sh
    swga sets -i primer_graph.gr -o found_sets.txt
    ```

6. The `swga score` command takes the output from `swga sets` and runs further statistics on the sets, only
    outputting a set if the primers in the set have a max distance on the foreground genome below a certain
    threshold. It then calculates a user-definable score for the set. The output of this command has the format
    per row of `[seq1,seq2,seq3...]    [score]     [{values used to compute score}]`.
    ```sh
    swga score -i found_sets.txt -o valid_sets.txt
    ```

There are often billions of valid sets possible, and it would be extremely time-consuming to find all of them. Therefore, `process_sets` can read the output from `find_sets` live and will exit after it finds `max_sets` valid primer sets. Alternatively, one could let the `find_sets` write to a file and find all sets, then pass it to `process_sets` and test different parameters.


## Theory

We need to select a compatible set of primers that binds to a target (foreground) genome with a higher frequency than a background genome for selective whole genome amplification. The primers need to be compatible with each other- i.e. not heterodimers- and have appropriate foreground and background binding frequencies.

Given a list of primers, the goal of finding a group of _n_ compatible primers (where _n_ is the desired set size) is mappable to finding the largest [clique](https://en.wikipedia.org/wiki/Clique_(graph_theory)) in a graph.

A clique in graph theory is a set of nodes in an undirected graph that are all connected to one another. If we visualize all the primers as nodes on a graph, and connect each compatible pair of primers with an edge, we can see that a set of mutually compatible primers takes the form of a clique in this graph.

Unfortunately, this is a known NP-complete problem, so no known efficient solution exists. Exponential-time solutions do exist, however, and we will be using an implementation of the branch-and-bound algorithm called [cliquer](http://users.tkk.fi/~pat/cliquer.html) to find the desired primer sets.

# Contributors and Acknowledgements

## Coding:
  - [Erik Clarke](https://github.com/eclarke) -- pipeline, filtering, set finding
  - [Calvin Morrison](https://github.com/mutantturkey) -- heterodimer comparison, k-mer counting
