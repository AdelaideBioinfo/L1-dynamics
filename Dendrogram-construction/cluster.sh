#!/bin/bash

# Cluster hits for species with lots of L1s (e.g. mammals), before making dendrograms

# Invoked by:
#
# GENOME=Tupaia.belangeri ORDER=Mammalia CLUSTERNO=7 sbatch cluster.sh
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1:00:00 
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au 

# mkdir -p UpperCase
# Convert all lower case to upper case (ignoring header)
# do this because usearch cluster_fast ignores masking parameters
#cat "$GENOME"_L1_annotated_"$CUTOFF".fasta \
#| awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' \
#> UpperCase/"$GENOME"_upper.fasta
#cd UpperCase

# go to dir
cd /data/rc003/atma/LASTZ_extraction/Annotated_L1s/$ORDER/$GENOME/UpperCase/active_clusters/

# Cluster the seqs at 80% identity
# Change the -id parameter to cluster at higher identites
usearch -cluster_fast "$GENOME"_cluster_"$CLUSTERNO" -id 0.80 --maxrejects 100 --maxaccepts 10 -uc "$GENOME"_clusters.uc -clusters "$GENOME"_cluster_"$CLUSTERNO"_
