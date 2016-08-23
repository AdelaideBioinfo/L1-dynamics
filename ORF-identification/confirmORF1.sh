#!/bin/bash

# HMMer confirmation for amino acid ORF1 candidates

# $1 = species name
# $2 = order
# Example usage: ./confirmORF1.sh Homo.sapiens Mammalia

###### Move to the desired directory ######
# since there are a lot of output files, 
# put each species in a separate directory
cd /data01/Protein_dbs/ORF_validation/$2/$1

###### Filter unknown characters before scanning ######
# replace '?' with 'X'
sed -i 's/?/X/g' "$1"_extended_ORF1_candidates_translation.fasta

###### Scan candidates (queries) against known protein domains ######
# using the Pfam-A.hmm database (note: this needs to be downloaded from Pfam website)
hmmscan --domtblout "$1"_extended_ORF1_candidates_translation.fasta.out /data01/Protein_dbs/Pfam-A.hmm \
"$1"_extended_ORF1_candidates_translation.fasta > "$1"_extended_ORF1_candidates_translation.fasta.log

###### Make a list of all the query seq names (for reference) ######
# list all the headers from the orf candidates that were queried against Pfam
grep '>' "$1"_extended_ORF1_candidates_translation.fasta | sed 's/>//g' > "$1"_ORF1_queries.txt
# remove empty lines from this file
sed -i '/^$/d' "$1"_ORF1_queries.txt

###### Discard ORF1 candidates that are actually ORF2 ######
# done by recognising RT and Endo domains
# first, remove comment lines
sed -i '/^#/ d' "$1"_extended_ORF1_candidates_translation.fasta.out
# extract the fields of interest
awk '{print $1 "\t" $4}' "$1"_extended_ORF1_candidates_translation.fasta.out > "$1"_ORF1_filtered.tmp
# separate each ORF candidate sequence by a blank line
# ie. group all the domains found in each candidate seq
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' "$1"_ORF1_filtered.tmp > "$1"_ORF1_sep.tmp
# grab the headers from output that contain rt or endo domains
# these are the ones we wish to exclude
cat "$1"_ORF1_sep.tmp | cut -f 1,2 \
| awk '{if ( $1=="RVT_1" || $1=="RVT_3" || $1=="Exo_endo_phos" || $1=="Exo_endo_phos_2" || $1=="zf-RVT" ) print $2 }' \
| sort -u > "$1"_actuallyORF2.tmp

###### Find the candidates that do not contain any ORF2 domains ######
# i.e. list headers excluding those in "$1"_actuallyORF2.tmp
comm -2 -3 <(sort "$1"_ORF1_queries.txt) <(sort "$1"_actuallyORF2.tmp) \
> "$1"_ORF1_candidates_excluding_ORF2.tmp
# this gives a much smaller, more likely ORF1 subset to screen for known domains

###### Extract FASTA sequences of the new subset of ORF1 candidates ######
# extract amino acid seqs from the original query FASTA file
cat "$1"_extended_ORF1_candidates_translation.fasta \
| awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' \
| awk -F"\t" 'BEGIN{while((getline k < "'$1'_ORF1_candidates_excluding_ORF2.tmp")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' \
> "$1"_ORF1_candidates_excluding_ORF2.fasta.tmp

###### Rescan this new subset against Pfam database ######
hmmscan --domtblout "$1"_ORF1_candidates_excluding_ORF2.fasta.out.tmp /data01/Protein_dbs/Pfam-A.hmm \
"$1"_ORF1_candidates_excluding_ORF2.fasta.tmp > "$1"_ORF1_candidates_excluding_ORF2.fasta.log.tmp

###### Confirm ORF1 candidates by looking for common domains ######
# remove comment lines
sed -i '/^#/ d' "$1"_ORF1_candidates_excluding_ORF2.fasta.out.tmp
# extract the fields of interest
awk '{print $1 "\t" $4}' "$1"_ORF1_candidates_excluding_ORF2.fasta.out.tmp > "$1"_ORF1_filtered2.tmp
# separate each ORF candidate sequence by a blank line
# ie. group all the domains found in each candidate seq
awk -v i=2 'NR>1 && $i!=p { print "" }{ p = $i } 1' "$1"_ORF1_filtered2.tmp > "$1"_ORF1_sep2.tmp
# grab the headers from output that contain known ORF1 domains
cat "$1"_ORF1_sep2.tmp | cut -f 1,2 \
| awk '{if ( $1=="Transposase_22" || $1=="zf-C2H2_10" || $1=="zf-C2H2_2" || $1=="RRM_1" || $1=="RRM_3" || $1=="RRM_5" || $1=="RRM_6" || $1=="RRM_occluded" || $1=="Nup35_RRM_2" || $1=="RRM_7" || $1=="zf-CCHC" || $1=="zf-CCHC_2" || $1=="zf-CCHC_3" || $1=="zf-CCHC_4" || $1=="zf-CCHC_5" || $1=="zf-CCHC_6" ) print $2 }' \
| sort -u > "$1"_confirmedORF1.txt

###### Extract FASTA sequences from confirmed ORF1 ######
# extract the amino acid seq of the confirmed ORFs from the original query FASTA file
cat "$1"_extended_ORF1_candidates_translation.fasta \
| awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' \
| awk -F"\t" 'BEGIN{while((getline k < "'$1'_confirmedORF1.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' \
> "$1"_confirmedORF1.fasta

###### Keep track of the candidates without confirmed ORF1 domains ######
# list headers that did not contain confirmed seqs
# (either had no domains, or no 'known' domains)
comm -2 -3 <(sort "$1"_ORF1_candidates_excluding_ORF2.tmp) <(sort "$1"_confirmedORF1.txt) > "$1"_unconfirmedORF1.txt

###### Keep track of the domains found in each ORF1 candidate sequence ######
# (both confirmed and unconfirmed - but not candidates which were actually ORF2)
# merge domains according to the header name
# so that each line has a unique header
awk '$2 != prev {if (NR != 1) print prev; prev=$2; delete a};
!($1 in a){a[$1]++; printf "%s,", $1};
END {print prev}' "$1"_ORF1_filtered2.tmp > "$1"_ORF1_merged.tmp
# rearrange to put into nicer format
# e.g. extract all domains and seperate them by commas
cat "$1"_ORF1_merged.tmp \
| awk -F "," '{$NF = ""; print $0}' \
| sed 's/ /,/g' \
| sed 's/\,$//g' > 1sthalf.tmp
# then extract the corresponding seq name
cat "$1"_ORF1_merged.tmp \
| awk -F "," '{print $NF}' > 2ndhalf.tmp
# paste together (tab-delimited) with col1=header, col2=domains
paste 2ndhalf.tmp 1sthalf.tmp > "$1"_ORF1_domains.txt

###### Write out a summary log file ######
# print the original number of ORF1 candidate seqs
wc -l "$1"_ORF1_queries.txt \
| awk '{print "# original ORF1 candidates: " $1 "\n"}' \
> "$1"_ORF1_summary.txt
# print how many turned out to be ORF2
wc -l "$1"_actuallyORF2.tmp \
| awk '{print "# that were actually ORF2: " $1 "\n"}' \
>> "$1"_ORF1_summary.txt
# print the new subset of ORF1 candidates
wc -l "$1"_ORF1_candidates_excluding_ORF2.tmp \
| awk '{print "# new subset of ORF1 candidates (excluding ORF2): " $1}' \
>> "$1"_ORF1_summary.txt
# print how many confirmed ORF1 there are
wc -l "$1"_confirmedORF1.txt \
| awk '{print "\n# candidates with confirmed domains: " $1}' \
>> "$1"_ORF1_summary.txt
# check that the FASTA file of confirmed ORF1 is the same
grep -c '>' "$1"_confirmedORF1.fasta \
| awk '{print "\n# extracted FASTA seqs (should be same): " $1}' \
>> "$1"_ORF1_summary.txt
# print the domains in the confirmed ORF1 seqs
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_confirmedORF1.txt "$1"_ORF1_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/\nORF1 domains:\n/' \
>> "$1"_ORF1_summary.txt
# print how many ORF1 do not contain known domains
wc -l "$1"_unconfirmedORF1.txt \
| awk '{print "\n# candidates that do not contain known ORF1 domains: " $1}' \
>> "$1"_ORF1_summary.txt
# print the domains in the unconfirmed ORF1 candidates
# not including RVT or Endo domains
awk 'NR==FNR{c[$1]++;next};c[$1]>0' "$1"_unconfirmedORF1.txt "$1"_ORF1_domains.txt \
| awk '{print $2}' \
| tr "," "\n" \
| grep -v "^\s*$" \
| sort \
| uniq -c \
| sort -bnr \
| awk '{print $1 " " $2}' \
| sed '1s/^/\nNon-ORF1 domains:\n/' \
>> "$1"_ORF1_summary.txt

###### Remove unnecessary files ######
# i.e. all files with extension .tmp
rm *.tmp

# Should be left with 9 files in the directory:
# 1. "$1"_extended_ORF1_candidates_translation.fasta: original query candidate seqs, amino acid form
# 2. "$1"_extended_ORF1_candidates_translation.fasta.out: raw hmmscan output
# 3. "$1"_extended_ORF1_candidates_translation.fasta.log: hmmscan log file
# 4. "$1"_ORF1_queries.txt: header names of the queries/candidates
# 5. "$1"_ORF1_domains.txt: tab-delimited file of the domains found in each ORF1 candidate
# note: this file only contains real ORF1 candidates (not those that turned out to be ORF2)
# 6. "$1"_confirmedORF1.txt: header names of all confirmed ORF1s
# 7. "$1"_confirmedORF1.fasta: FASTA seqs of all confirmed ORF1s
# 8. "$1"_unconfirmedORF1.txt: header names of all unconfirmed ORF1s 
# (unconfirmed from the real ORF1 candidates - does not include any that were found to be ORF2)
# 9. "$1"_ORF1_summary.txt: summary numbers showing 
				# how many original candidates in this species,
				# how many turned out to be ORF2,
				# the new subset of ORF1 candidates after removing ORF2 seqs,
                # how many confirmed ORF1 (containing known domains),
                # (and a sanity check of the FASTA file)
                # the domains within confirmed ORF1 seqs
                # how many unconfirmed ORF1 (not ORF2, but don't contain known domains either)
                # domains within the unconfirmed ORF1
