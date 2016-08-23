#!/bin/bash

# Used to extract the nucleotide sequences of the TBLASTN hits
# Can then use these nucleotide L1s as further queries with LASTZ
# To search the genomic data

# Invoked by:
#
# sbatch getNuclSeq.sh
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1:00:00 
#SBATCH --mem=30GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au      

# Load the necessary modules
module load BEDTools/2.25.0-foss-2015b

# Extract nucleotide FASTA sequence for each nt and htgs L1 hit
# Again, this example uses the two amphibian species taxid ids
arr=( 'txid125878' 'txid8364' )
for i in "${arr[@]}"
do

# Grab the headers from the amino acid FASTA file
grep '>' "$i"_sseq.fasta > "$i"_sseq_headers.txt

# Pull out all the plus strand sequences
cat "$i"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| sed 's/-/\t/g' \
| awk '{if ($2 < $3) print $1 "\t" $2 "\t" $3 "\t" "L1hit" "\t" "1" "\t" "+"}' \
> plus.txt 

# Pull out all the minus strand sequences
cat "$i"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| sed 's/-/\t/g' \
| awk '{if ($2 > $3) print $1 "\t" $3 "\t" $2 "\t" "L1hit" "\t" "1" "\t" "-"}' \
> minus.txt

# Combine minus and plus strand sequences
cat minus.txt plus.txt > all.txt

# Sort by chr/scaffold and then by start position in ascending order 
bedtools sort -i all.txt > sorted.txt

# Merge nested or overlapping intervals
bedtools merge -s -i sorted.txt -c 4 -o collapse > merged.txt

# Merging broke the BED-like format
# (i.e. strand is now in 4th column, merged names in 5th column, etc)
# Rearrange the columns as required
awk '{print $1 "\t" $2 "\t" $3 "\t" "L1hit" "\t" "1" "\t" $4}' merged.txt > merged_"$i".bed

# Extract FASTA from corrected merged BED file
# Note: for this step, you will need to download the NCBI nt/htgs databases
bedtools getfasta -s -fi /data/rc003/atma/NCBI_databases/Both/nt_and_htgs.fasta -bed merged_"$i".bed -fo results.fasta

# Sort sequences by length
usearch -sortbylength results.fasta -fastaout "$i"_nucl_seqs.fasta

# Move to the allocated directory
mv "$i"_nucl_seqs.fasta ../Nucleotide_sequences/

# Remove unnecessary files
rm minus.txt
rm plus.txt
rm all.txt
rm sorted.txt
rm merged.txt
rm results.fasta

done
