#!/bin/bash

# Used to perform an all-against-all BLAST and Silix clustering 
# (e.g. for L1 ORF1 protein sequences)

# Make a FASTA file of all the hits of interest 
# Make sure every sequence has a unique name and includes the species name
# Format the FASTA file of seqs as a database 
formatdb -i L1_ORF1p_hits.fasta -p T -o T -n ORF1p
# -p: set to F for nucleotide seqs
# -o: set to T to parse seqID and create indexes
# -n: basename for this formatted database, e.g. BovBHits

# Perform an all-against-all blast of the hits in the database
blastall -p blastp -d ORF1p -i L1_ORF1p_hits.fasta -m 8 -e 1e-10 -a 8 -o blastall_ORF1p.out
# -p: chooses the BLAST program
# -d: database to query
# -i: query seqs (since we want to do an all-by-all comparison, the query seqs are the same as the database seqs)
# -m: for tabular output (to use with SiLiX), set to 8
# -e: e-value
# -a: number of CPUs (e.g. set to 8)
# -o: output file

# Use SiLiX software to cluster the hits (with default parameters)
silix L1_ORF1p_hits.fasta blastall_ORF1p.out > ORF1p.fnodes

# Alternatively, cluster the hits at a range of different % ids
for j in 60 65 70 75 80 85 90 95 98
do
silix L1_ORF1p_hits.fasta blastall_ORF1p.out -i "0.$j" -r 0.70 > ORF1p_$j.fnodes
done
# -i: min % identity
# -r: min length/overlap

# Sort the lines in each file numerically (e.g. family 1 to whatever-largest-number-is)
for j in 60 65 70 75 80 85 90 95 98
do
sort -n ORF1p_$j.fnodes > ORF1p_sorted_$j.fnodes
done

# Separate each family/cluster by a blank line
for j in 60 65 70 75 80 85 90 95 98
do
awk -v i=1 'NR>1 && $i!=p { print "" }{ p = $i } 1' ORF1p_sorted_$j.fnodes > ORF1p_separated_$j.fnodes
done
# i.e. on any line after the first, if the value of the "i"-th column (where i=1) is different to the previous value, print a blank line. 
# 1 at the end means true - so that awk prints the line

# Replace all '_' with spaces
for j in 60 65 70 75 80 85 90 95 98
do
sed "s/_/ /g" ORF1p_separated_$j.fnodes > ORF1p_columns_$j.fnodes
done

# Print all the different combinations of species found (ignoring species-specific clusters)
for j in 60 65 70 75 80 85 90 95 98
do
awk -f tst.awk ORF1p_columns_$j.fnodes > ORF1p_species_combinations_$j.fnodes
done
# output format is of the form: ClusterNo, Species
# Note: the ClusterNo shown is the FIRST cluster at which this unique species combination appears
# There may be later clusters which also show the same animal combination
# So if the animal combination is interesting, look for similar clusters after this ClusterNo
