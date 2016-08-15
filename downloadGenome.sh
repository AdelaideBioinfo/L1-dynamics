#!/bin/bash
# Shows how to download a publicly available genome from NCBI Assembly.
# To download multiple genomes at once, 
# create a for loop with the changing parameters being:
# directory name (here given as the species name, e.g. Elephantulus.edwardii)
# genome name (e.g. EleEdw1.0.fa)
# and ftp link (see below)

# Move to allocated directory
cd Genomes/

# Mammals
# scaffolded genome - Cape elephant shrew
mkdir -p Elephantulus.edwardii
cd Elephantulus.edwardii
wget --timestamping 'ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Elephantulus_edwardii/EleEdw1.0/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fa.gz'
gunzip unplaced.scaf.fa.gz
mv unplaced.scaf.fa EleEdw1.0.fa
cd ..
