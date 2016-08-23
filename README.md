# L1-dynamics
Sample code to accompany the L1 evolutionary dynamics across eukaryotes manuscript. Shows how to perform two independent extraction methods: iterative search using LASTZ on genomic data, versus translated nucleotide search of NCBI databases using TBLASTN. Subsequent analyses use programs such as MUSCLE, USEARCH, HMMER, etc.  

Order of execution:
L1-extraction (LASTZ, TBLASTN) -> ORF-identification -> Dendrogram-construction -> RT-identification -> Clustering-analysis

Within L1-extraction: 
#### LASTZ ####
downloadGenome.sh -> bundle.go -> renameToSeq.sh -> lastzExtractFromGenome.sh -> confirmLastzHits.sh
#### TBLASTN ####
tblastnExtractFromDatabase.sh -> getNuclSeq.sh -> confirmTblastnHits.sh -> rerun LASTZ pipeline

Within ORF-identi
