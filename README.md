PrimerSets
==========

Module to select sets of primers

Theory
----------

We need to select a compatible set of primers that binds to a target (foreground) genome with a higher frequency than a background genome for selective whole genome amplification. The primers need to be compatible with each other- i.e. not heterodimers- and have appropriate foreground and background binding frequencies. 

Given a list of primers, the goal of finding a group of _n_ compatible primers (where _n_ is the desired set size) is mappable to finding the largest [clique](https://en.wikipedia.org/wiki/Clique_(graph_theory)) in a graph.

A clique in graph theory is a set of nodes in an undirected graph that are all connected to one another. If we visualize all the primers as nodes on a graph, and connect each compatible pair of primers with an edge, we can see that a set of mutually compatible primers takes the form of a clique in this graph.

Unfortunately, this is a known NP-complete problem, so no known efficient solution exists. Exponential-time solutions do exist, however, and we will be using an implementation of the branch-and-bound algorithm called [cliquer](http://users.tkk.fi/~pat/cliquer.html) to find the desired primer sets.

Development will be focused on wrapping cliquer in Python to manage the workflow, devising a way of splitting the computational workload amongst multiple cores, and integrating with the rest of the SWGA pipeline.
