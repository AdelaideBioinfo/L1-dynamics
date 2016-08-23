#!/bin/bash

# Infer a maximum likelihood phylogeny from each alignment using FastTree

# Go to the directory containing the FASTA files
cd /home/Species_trees/Plants/

# Infer a maximum likelihood phylogeny using the alignment
# This tree can be visualised using tree viewer programs
fasttree -nt -gtr < L1_"$i"_seqs.afa > L1_"$i"_seqs.newick
