#!/bin/bash

# Invoked by:
#
# GENOME=Tarsius.syrichta ORDER=Mammalia CLUSTERNO=7 sbatch alignActiveClusters.sh
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1-00:00 
#SBATCH --mem=60GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au 

# load module muscle
module load MUSCLE/3.8.31

# go to dir
cd /data/rc003/atma/LASTZ_extraction/Annotated_L1s/$ORDER/$GENOME/UpperCase/active_clusters/

# align seqs in each 'active' cluster (i.e. contains at least one ORF2-intact seq)
muscle -in "$GENOME"_cluster_"$CLUSTERNO" -out "$GENOME"_cluster_"$CLUSTERNO".afa -maxiters 2 \
2> "$GENOME"_cluster_"$CLUSTERNO".log
