#!/bin/bash

# HMMer confirmation for amino acid ORF2 candidates

# $1 = species name
# $2 = order
# Example usage: ./confirmORF2.sh Homo.sapiens Mammalia

###### Move to the desired directory ######
# since there are a lot of output files, 
# put each species in a separate directory
cd /data01/Protein_dbs/ORF_validation/$2/$1

###### Filter unknown characters before scanning ######
# replace '?' with 'X'
sed -i 's/?/X/g' "$1"_extended_ORF2_candidates_translation.fasta

###### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database
hmmscan --domtblout "$1"_extended_ORF2_candidates_translation.fasta.out /data01/Protein_dbs/Pfam-A.hmm \
"$1"_extended_ORF2_candidates_translation.fasta > "$1"_extended_ORF2_candidates_translation.fasta.log

###### Make a list of all the query seq names (for reference) ######
# list all the headers from the orf candidates that were queried against Pfam
grep '>' "$1"_extended_ORF2_candidates_translation.fasta | sed 's/>//g' > "$1"_ORF2_queries.txt
# remove empty lines from this file
sed -i '/^$/d' "$1"_ORF2_queries.txt

###### Confirm ORF2 candidates by recognising RT domains ######
# remove comment lines
sed -i '/^#/ d' "$1"_extended_ORF2_candidates_translation.fasta.out
# extract the fields of interest
awk '{print $1 "\t" $4}' "$1"_extended_ORF2_candidates_translation.fasta.out > "$1"_ORF2_filtered.tmp
# separate each ORF candidate sequence by a blank line
# ie. group all the domains found in each candidate seq
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' "$1"_ORF2_filtered.tmp > "$1"_ORF2_sep.tmp
# grab the headers from output that contain "RVT_*" (ORF2 reverse transcriptase)
cat "$1"_ORF2_sep.tmp | cut -f 1,2 \
| awk '{if ( $1=="RVT_1" || $1=="RVT_3" ) print $2 }' \
| sort -u > "$1"_confirmedORF2.txt

###### Extract FASTA sequences from confirmed ORF2 ######
# extract the amino acid seq of the confirmed ORFs from the original query FASTA file
cat "$1"_extended_ORF2_candidates_translation.fasta \
| awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' \
| awk -F"\t" 'BEGIN{while((getline k < "'$1'_confirmedORF2.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' \
> "$1"_confirmedORF2.fasta

###### Keep track of the candidates without confirmed RT domains ######
# list headers that did not contain confirmed seqs
# (either had no recognisable domains, or none that contained RVT)
comm -2 -3 <(sort "$1"_ORF2_queries.txt) <(sort "$1"_confirmedORF2.txt) > "$1"_unconfirmedORF2.txt

###### Keep track of the domains found in each ORF2 candidate sequence ######
# (both confirmed and unconfirmed)
# merge domains according to the header name
# so that each line has a unique header
awk '$2 != prev {if (NR != 1) print prev; prev=$2; delete a};
!($1 in a){a[$1]++; printf "%s,", $1};
END {print prev}' "$1"_ORF2_filtered.tmp > "$1"_ORF2_merged.tmp
# rearrange to put into nicer format
# e.g. extract all domains and seperate them by commas
cat "$1"_ORF2_merged.tmp \
| awk -F "," '{$NF = ""; print $0}' \
| sed 's/ /,/g' \
| sed 's/\,$//g' > 1sthalf.tmp
# then extract the corresponding seq name
cat "$1"_ORF2_merged.tmp \
| awk -F "," '{print $NF}' > 2ndhalf.tmp
# paste together (tab-delimited) with col1=header, col2=domains
paste 2ndhalf.tmp 1sthalf.tmp > "$1"_ORF2_domains.txt

###### Write out a summary log file ######
# print the original number of ORF2 candidate seqs
wc -l "$1"_ORF2_queries.txt \
| awk '{print "# ORF2 candidates: " $1}' \
> "$1"_ORF2_summary.txt
# print how many confirmed ORF2 there are
wc -l "$1"_confirmedORF2.txt \
| awk '{print "\n# candidates that contain RVT: " $1}' \
>> "$1"_ORF2_summary.txt
# check that the FASTA file of confirmed ORF2 is the same
grep -c '>' "$1"_confirmedORF2.fasta \
| awk '{print "\n# extracted FASTA seqs (should be same): " $1}' \
>> "$1"_ORF2_summary.txt
# print the domains in the confirmed ORF2 seqs
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_confirmedORF2.txt "$1"_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/\nORF2 domains:\n/' \
>> "$1"_ORF2_summary.txt
# print how many ORF2 do not contain RVT
wc -l "$1"_unconfirmedORF2.txt \
| awk '{print "\n# candidates that do not contain RVT: " $1}' \
>> "$1"_ORF2_summary.txt
# print the domains in the unconfirmed ORF2 candidates
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_unconfirmedORF2.txt "$1"_ORF2_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/non-ORF2 domains:\n/' \
>> "$1"_ORF2_summary.txt


###### Remove unnecessary files ######
# i.e. all files with extension .tmp
rm *.tmp

# Should be left with 9 files in the directory:
# 1. "$1"_extended_ORF2_candidates_translation.fasta: original query candidate seqs, amino acid form
# 2. "$1"_extended_ORF2_candidates_translation.fasta.out: raw hmmscan output
# 3. "$1"_extended_ORF2_candidates_translation.fasta.log: hmmscan log file
# 4. "$1"_ORF2_queries.txt: header names of the queries/candidates
# 5. "$1"_ORF2_domains.txt: tab-delimited file of the domains found in each ORF candidate
# 6. "$1"_confirmedORF2.txt: header names of all confirmed ORF2s
# 7. "$1"_confirmedORF2.fasta: FASTA seqs of all confirmed ORF2s
# 8. "$1"_unconfirmedORF2.txt: header names of all unconfirmed ORF2s
# 9. "$1"_ORF2_summary.txt: summary numbers showing
							# how many original candidates in this species,
                            # how many confirmed ORF2 (containing RT),
                            # (and a sanity check of the FASTA file)
                            # domains within confirmed ORF2 seqs
                            # how many unconfirmed ORF2
                            # and the domains within unconfirmed ORF2
