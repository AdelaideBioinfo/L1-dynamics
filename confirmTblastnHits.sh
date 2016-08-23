#!/bin/bash

# Confirms which TBLASTN hits are really L1 sequences
# Using a reciprocal best hit check
# All extracted nucleotide sequences are screened using CENSOR,
# against the Repbase library of known repeats
# Hits are kept if the best hit from CENSOR is an L1
# and discarded if the best hit is another repeat (e.g. CR1, hAT)

# Invoked by:
#
# sbatch confirmTblastnHits.sh
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1-00:00 
#SBATCH --mem=30GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au      

# Load the necessary modules
module load bioperl/1.6.924
module load censor/4.2.29
module load wu-blast/2.0
module load BEDTools/2.25.0-foss-2015b

# This array uses the two amphibian species as examples
arr=( 'txid125878' 'txid8364' )
for i in "${arr[@]}"
do

# Run CENSOR on the extracted nucleotide seqs, using the Repbase library
# Note: download the latest Repbase library from http://www.girinst.org/repbase/
censor -bprm cpus=8 -lib /data/rc003/atma/tblastn/Nucleotide_sequences/RepBase21.03_all_seqs.ref "$i"_nucl_seqs.fasta

# Want to keep only the nucleotide sequences that are verified as L1s

# First, extract map file lines in the right orientation (d, not c)
# Print the header (which contains coords and strand) and column 4 (hit name)
# Put into BED-like format
cat "$i"_nucl_seqs.fasta.map \
| awk '{if ($7=="d") print $0 }' \
| awk '{print $1 "\t" $4}' \
| sed 's/:/\t/g' \
| sed 's/(-)/.rev/g' \
| sed 's/(+)/.fwd/g' \
| awk '{gsub(/-/,"\t",$2); print}' \
| sed 's/.rev/\t-/g' \
| sed 's/.fwd/\t+/g' \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' \
> rearranged.map

# Sort the file by chr/scaffold then coordinates
bedtools sort -i rearranged.map > sorted.map

# Merge overlapping hits
bedtools merge -s -i sorted.map -c 4 -o collapse > merged.map

# Rearrange columns 
awk '{print $5 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' merged.map > merged.txt

# Extract sequences which identified as L1s
# Download from Repbase all potential L1 sequences
# e.g. everything in L1 group and Tx1 group (1826 seqs total)
# Extract all the header names -> 1826_L1_and_Tx1_repbase_headers.txt
# Compare this to the hit names in the merged.txt file
# and then rearrange into BED format again
awk 'NR==FNR { a[$1] = $1; next} { for (k in a) if ($1 ~ a[k]) { print $0; break } }' /data/rc003/atma/tblastn/Nucleotide_sequences/1826_L1_and_Tx1_repbase_headers.txt merged.txt \
| awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5}' \
> "$i"_L1.bed

# Extract FASTA from L1-only merged BED file
bedtools getfasta -s -fi /data/rc003/atma/NCBI_databases/Both/nt_and_htgs.fasta -bed "$i"_L1.bed -fo results.fasta

# Sort sequences by length
usearch -sortbylength results.fasta -fastaout "$i"_L1.fasta

# Remove unnecessary files
rm rearranged.map
rm sorted.map
rm merged.map
rm merged.txt
rm results.fasta

done
