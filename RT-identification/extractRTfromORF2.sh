#!/bin/bash

# Extract RT domain from confirmed ORF2 sequences

# $1 = species name
# $2 = order
# Example usage: ./extractRTfromORF2.sh Homo.sapiens Mammalia

###### Move to the desired directory ######
# since there are a lot of output files, 
# put each species in a separate directory
cd /data01/Protein_dbs/ORF_validation/$2/$1

###### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database
hmmscan --domtblout "$1"_confirmedORF2.fasta.out /data01/Protein_dbs/Pfam-A.hmm \
"$1"_confirmedORF2.fasta \
> "$1"_confirmedORF2.fasta.log

# pull out a BED file containing the location of the RT domains
# check that the RT domains validates as RVT_1 or RVT_3
# and that it is full-size (e.g. expect it to be around 255 amino acids,
# so set a restriction of >200 amino acids
# note: $20 and $21 are the envelope coordinates on the query protein seq
# so even if the alignment doesn't stretch as far, or scores low in some regions
# this 'envelope' spans the entire suspected RT domain
cat "$1"_confirmedORF2.fasta.out \
| awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' \
> "$1"_RT_domain.bed

# retrieve the fasta using this BED file
# check that you are not off by 1!! (e.g. if HMMer is 1-based and BEDTools is 0-based)
# no need to use strand since they are all facing the same way
fastaFromBed -fi "$1"_confirmedORF2.fasta -bed "$1"_RT_domain.bed -fo "$1"_RT_domain.fasta

# move to a easy-to-find location
mv "$1"_RT_domain.fasta /data01/Protein_dbs/RT_domains/
