#!/bin/bash

# Used to extract L1 sequences from each genome
# Requires as input the query L1 file and genome chromosomes/scaffolds
# Very large genomes need to be divided into smaller files for input to LASTZ (e.g. use bundle.go)

# $1 = order (e.g. Ecdysozoa)
# $2 = genome directory name (e.g. Anopheles.gambiae)
# $3 = startvalue (e.g. 1)
# $4 = endvalue (e.g. 7)
# $5 = query file (e.g. 1826_L1_and_Tx1_repbase.fasta)

# Usage: 
# lastzExtractFromGenome.sh Ecdysozoa Anopheles.gambiae 1 7 1826_L1_and_Tx1_repbase.fasta

cd /data/rc003/atma/LASTZ_extraction/Results
# Make directory for this species order
mkdir -p $1
cd $1
# Make directory for this specific genome
mkdir -p $2
# Make sub-directory to store L1/Order/Genome alignment hits
cd $2
mkdir -p Hits
cd Hits

# Align L1 query seqs to all Genome seqs, using LASTZ
# Note that this currently points to the Split_seqs genome sub-directory
seq $3 $4 | parallel lastz '/data/rc003/atma/LASTZ_extraction/Genomes/'$1'/'$2'/Split_seqs/seq{}.fa[unmask,multiple] /data/rc003/atma/LASTZ_extraction/Query/'$5'[unmask,multiple] --chain --gapped --coverage=80 --ambiguous=n --ambiguous=iupac --format=general-:name2,start2,end2,score,strand2,size2,name1,start1,end1 > /data/rc003/atma/LASTZ_extraction/Results/'$1'/'$2'/Hits/LASTZ_L1_'$2'_'$5'_seq{}'

# Make sure you are in the right directory
cd /data/rc003/atma/LASTZ_extraction/Results/$1/$2/Hits

# Remove all files that are empty
find -size  0 -print0 | xargs -0 rm

# Concatenate all hit files into one file
cat LASTZ_L1_"$2"_"$5"_seq* > LASTZ_L1_"$2"_"$5"_AllSeqs

# Rearrange columns to put the concatenated file in BED-like form (for BEDTools) 
awk '{print $7 "\t" $8 "\t" $9 "\t" $1 "\t" "1" "\t" $5}' LASTZ_L1_"$2"_"$5"_AllSeqs >> BedFormat_L1_"$2"_"$5"_AllSeqs

# Sort by chr/scaffold and then by start position in ascending order 
bedtools sort -i BedFormat_L1_"$2"_"$5"_AllSeqs > Sorted_L1_"$2"_"$5"_AllSeqs

# Merge nested or overlapping intervals
bedtools merge -s -i Sorted_L1_"$2"_"$5"_AllSeqs -c 4 -o collapse > Merged_L1_"$2"_"$5"_AllSeqs

# Merging broke the BED-like format
# (i.e. strand is now in 4th column, merged names in 5th column, etc)
# Rearrange the columns as required
awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' Merged_L1_"$2"_"$5"_AllSeqs >> Merged_L1_"$2"_"$5"_AllSeqs_proper.bed

# Extract FASTA from corrected merged BED file
# Note that (for simplicity) this uses the whole genome file, not the Split_seqs
bedtools getfasta -s -fi /data/rc003/atma/LASTZ_extraction/Genomes/$1/$2/*.fa -bed Merged_L1_"$2"_"$5"_AllSeqs_proper.bed -fo results.fasta

# Sort sequences by length
usearch -sortbylength results.fasta -fastaout "$2"_L1_"$5"

# Remove unnecessary files
rm results.fasta

# Move this final FASTA file to a separate directory (for easy access)
mkdir -p /data/rc003/atma/LASTZ_extraction/Results/$1/$2/FASTA
mv "$2"_L1_"$5" ../FASTA
