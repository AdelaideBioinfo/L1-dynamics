#!/bin/bash

# Invoked by:
#
# QUERY=/data/rc003/atma/tblastn/L1_query_ORFps.fasta sbatch tblastn.sh
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
module load BLAST+/2.2.31-foss-2015b-Python-2.7.11

# Run a remote tblastn (see tblastn -help for details)
# Specify the organism of interest using -entrez_query
# To do this, make an array of the taxid ids to be searched
# e.g. this is set to the two amphibian species, Nanorana parkeri and Xenopus tropicalis
arr=( 'txid125878' 'txid8364' )
for i in "${arr[@]}"
do

tblastn -db "nt htgs" -remote -entrez_query "$i"[Organism] -query "$QUERY" -out ../tblastn/Amphibia/"$i".out -outfmt "6 qseqid sseqid bitscore qcovs evalue pident length mismatch gapopen qstart qend qseq sstart send sseq" -evalue 1e-5 -show_gis \
2> ../tblastn/Amphibia/"$i".log
	
cat ../tblastn/Amphibia/"$i".out \
| sort -nr -k3 \
> ../tblastn/Amphibia/"$i"_sorted.out

# Put the target seq hits in FASTA format (for further checking)
# first, remove gaps ("-") from the target sequence (column 15)
# then print out all of the target seqs, with the header in the form: 
# >sseqid:sstart-send	qseqid:qstart-qend	score:*,qcov:*,evalue:*,pident:*
# e.g. >gi|699251403|ref|XM_002119377.3|:618-821	L1-1a_Cis_ORFp:208-291	score:38.5,qcov:17,evalue:6e-56,pid:30.68
cat ../tblastn/Amphibia/"$i"_sorted.out \
| awk -F" " '{gsub(/[-]/,"",$15)}1' OFS=" " \
| awk '{print ">" $2 ":" $13 "-" $14 "\t" $1 ":" $10 "-" $11 "\t" "score:" $3 ",qcov:" $4 ",eval:" $5 ",pid:" $6 "\n" $15}' \
> ../tblastn/Amphibia/"$i"_sseq.fasta
# Note: "$i"_sseq.fasta contains the protein hit sequences, not the nucleotide sequences

done
