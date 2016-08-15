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
