# Selective whole genome amplification
[![Build Status](https://travis-ci.org/eclarke/swga.svg?branch=master)](https://travis-ci.org/eclarke/swga)
## Introduction 
This is an easy-to-use, start-to-end package for finding sets of primers that selectively amplify a particular genome (the "foreground" genome) over a background genome. For instance, we can design a set of primers that amplify a parasite's genome in a sample that is overwhelmingly composed of host DNA.

You can run SWGA on hardware ranging from a Mac laptop to a high-end server. 

## Features:
- Counts all the possible primers in a size range in both genomes
- Filters primers based on:
  - foreground and background genome binding rates
  - melting temperatures (with a built-in melt temp calculator that accounts for mono- and divalent cation solutions!)
  - Possible homodimerization
- Finds primer sets containing primers that are compatible with each other using graph theory (largest clique formation). The process ensures:
  - No primer in a set is a heterodimer
  - Even binding site spacing in foreground genome
  - Low total binding to background genome
- Score each set based on certain binding metrics and allows exploration of high-scoring sets via output to common formats.

## Installation
Follow the installation instructions [here](https://github.com/eclarke/swga/wiki/Installation)

## Using SWGA
Follow the guide on our Wiki/[Quick Start](https://github.com/eclarke/swga/wiki) to get started!

## Updates
New features and bugfixes are released all the time. To update, simply follow steps 3-5 on the [installation instructions](https://github.com/eclarke/swga/wiki/Installation).
