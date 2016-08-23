#!/bin/bash

# Confirms L1 sequences found in the genome
# Using CENSOR to perform a recirprocal best hit check, as before

# Invoked by:
#
# ORDER=Ecdysozoa GENOME=Anopheles.gambiae sbatch confirmLastzHits.sh
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH --time=3-00:00 
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au      

# Load the necessary modules
module load bioperl/1.6.924
module load censor/4.2.29
module load wu-blast/2.0
module load BEDTools/2.25.0-foss-2015b

# Go to the right directory
cd /data/rc003/atma/LASTZ_extraction/testResults/$ORDER/$GENOME/FINAL/

# Screen the L1 candidates in this genome against the entire Repbase library of repeats
censor -bprm cpus=16 \
-lib /data/rc003/atma/tblastn/Nucleotide_sequences/RepBase21.03_all_seqs.ref \
/data/rc003/atma/LASTZ_extraction/testResults/$ORDER/$GENOME/FINAL/"$GENOME"_L1_candidates.fasta
# Outputs a *.map file in the FINAL directory

# Want to keep only the nucleotide sequences that are verified as L1s, as before
# So extract map file lines in the right orientation (d, not c)
# Print the header (which contains coords and strand) and column 4 (hit name)
# Put into BED-like format
cat "$GENOME"_L1_candidates.fasta.map \
| awk '{if ($7=="d") print $0}' | awk '{print $1 "\t" $4}' \
| sed 's/:/\t/g' \
| sed 's/(-)/.rev/g' | sed 's/(+)/.fwd/g' \
| awk '{gsub(/-/,"\t",$2); print}' \
| sed 's/.rev/\t-/g' | sed 's/.fwd/\t+/g' \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' \
> rearranged.map

# Sort the file by chr/scaffold then coordinates
bedtools sort -i rearranged.map > sorted.map

# Merge overlapping hits
bedtools merge -s -i sorted.map -c 4 -o collapse > merged.map

# Rearrange columns 
awk '{print $5 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' merged.map > merged.txt

# Check the number of hits after merging
wc -l merged.txt | awk '{print "merged.txt: " $1}'

# Extract sequences which identify as L1s
awk 'NR==FNR { a[$1] = $1; next} { for (k in a) if ($1 ~ a[k]) { print $0; break } }' \
/data/rc003/atma/tblastn/Nucleotide_sequences/1826_L1_and_Tx1_repbase_headers.txt merged.txt \
| awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5}' \
> "$GENOME"_L1_verified.bed

# Check the number of L1-verified hits
wc -l "$GENOME"_L1_verified.bed | awk '{print "L1_verified.bed: " $1}'

# Extract FASTA from L1-only merged BED file
bedtools getfasta -s -fi /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa \
-bed "$GENOME"_L1_verified.bed -fo output.fasta

# Sort sequences by length
usearch -sortbylength output.fasta -fastaout "$GENOME"_L1_verified.fasta

# Check the number of L1 sequences (should be the same as in L1-verified bed file)
grep -c '>' "$GENOME"_L1_verified.fasta | awk '{print "check FASTA file (should be same): " $1}'

# Move to a new directory
mkdir -p ../VERIFIED
mv "$GENOME"_L1_verified.fasta ../VERIFIED

# Remove unnecessary files
rm rearranged.map sorted.map merged.map merged.txt output.fasta
