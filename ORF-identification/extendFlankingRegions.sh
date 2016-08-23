#!/bin/bash

# Extend all L1 sequences by 1 kb (each side), 
# to analyse conservation of open reading frames

# Invoked by:
#
# ORDER=Ecdysozoa GENOME=Anopheles.gambiae sbatch extendFlankingRegions.sh
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

# Load modules
module load BEDTools/2.25.0-foss-2015b

# Go to directory
cd /data/rc003/atma/LASTZ_extraction/testResults/$ORDER/$GENOME/FINAL/

# Use the merged BED file from the extraction pipeline
# to extend each sequence by 1 kb either side 
bedtools slop -i "$GENOME"_L1_verified.bed -g /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa.fai -b 1000 > "$GENOME"_extended.bed

# Extract a FASTA file of the extended hits
bedtools getfasta -s -fi /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa -bed "$GENOME"_extended.bed -fo ext.fasta

# Sort the extended results by sequence length
usearch -sortbylength ext.fasta -fastaout "$GENOME"_extended.fasta

# Move FASTA file of extended hits to verified directory
mv "$GENOME"_extended.fasta ../VERIFIED

# Make a new directory
mkdir -p ../Tracking_names

# Copy important BED files
cp "$GENOME"_L1_verified.bed ../Tracking_names
cp "$GENOME"_extended.bed ../Tracking_names

# Move to the new directory
cd ../Tracking_names

# Change both original and extended BED files from bed format to hit format
cat "$GENOME"_L1_verified.bed | awk '{print $1 ":" $2 "-" $3 "(" $6 ")"}' > "$GENOME"_original_headers.txt
cat "$GENOME"_extended.bed | awk '{print $1 ":" $2 "-" $3 "(" $6 ")"}' > "$GENOME"_extended_headers.txt

# Paste original and extended files together;
# to keep a record of the changes, for annotation during dendrogram analysis
paste "$GENOME"_original_headers.txt "$GENOME"_extended_headers.txt > "$GENOME"_original_and_extended_headers.txt
#sed -i '1s/^/#original\t#extended\n/' "$GENOME"_original_and_extended_headers.txt
