# L1-dynamics
Sample code to accompany the L1 evolutionary dynamics across eukaryotes manuscript. Shows how to perform two independent extraction methods: iterative search using LASTZ on genomic data, versus translated nucleotide search of NCBI databases using TBLASTN. Subsequent analyses use programs such as MUSCLE, USEARCH, HMMER, etc.  

Order of execution:
L1-extraction (LASTZ, TBLASTN) -> ORF-identification -> Dendrogram-construction
                                                     -> RT-identification
                                                     -> Clustering-analysis
