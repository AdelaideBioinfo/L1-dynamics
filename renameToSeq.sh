#!/bin/bash

# Used to rename the output from bundle.go into appropriate format for lastz.sh

# Parameters:
# $1 = GENOME file that has been split using bundle.go

# Example usage:
# Put this script in ~/bin and chmod ugo+rwx renameToSeq.sh
# Run: renameToSeq.sh genome.fa

# rename genome.fa-0.fa to seq1.fa, genome.fa-1.fa to seq2.fa, etc
# first, ls all the genome.fa-*.fa files
# use "-" as a delimiter to get *.fa only (e.g. 0.fa, 1.fa, 2.fa, ...)
# sort numerically
# print the original file name in column 1, and the new file name in column 2
# rename everything in column 1 to the respective field in column 2
ls *.fa-*.fa | cut -d "-" -f 2 | sort -g | awk '{print "'$1'" "-" $0, "seq" NR ".fa"}' | xargs -n2 mv

# move seq* files to a new sub-directory
mkdir -p Split_seqs
mv seq*.fa Split_seqs

# print the sequence files in numerical order
# e.g. to see how many files need to be analysed
cd Split_seqs
ls seq* | cut -c 4- | sort -g | awk '{print "seq" $0}'
