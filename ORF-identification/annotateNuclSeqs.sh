#!/bin/bash

# Add ORF annotation to nucleotide FASTA L1 seqs (un-extended, original hits)

# $1 = species name
# $2 = order
# Example usage: annotateNuclSeqs.sh Homo.sapiens Mammalia

# Note: this was run on Backup harddrive. Change path directories to run on VM.

###### Move to the desired directory ######
# since there are a lot of output files, 
# put each species in a separate directory
cd /Volumes/Backup/ORF_validation/$2/$1

# grab the header names of the confirmed ORF1 seqs
# and annotate them as having ORF1
# note: these are the extended headers
# in case this file is blank, add a #comment
# and remove empty lines
cat "$1"_confirmedORF1.txt \
| awk -F "_-_ORF" '{print $1 "\t" "ORF1_"}' \
| awk 'BEGIN{print "#comment"}{print}END{print ""}' \
| sed '/^$/d' \
> "$1"_confirmedORF1_extendedheaders.tmp  

# already have the original_and_extended record from extend_flanking_regions.sh
# need to append ORF1 or ORF2 labels to this file
#awk 'FNR==NR {a[$1]=$2;next} {print $0,a[$2]}' \
awk 'FNR==NR {a[$1]=$2;next} {if ($2 in a) print $0,a[$2]; else print $0;}' \
"$1"_confirmedORF1_extendedheaders.tmp /Volumes/Backup/Results/$2/Tracking_names/"$1"_original_and_extended_headers.txt \
> "$1"_annotated_with_ORF1.tmp
# with awk comparison of two files, first {} refers to file1; second {} to file2
# a[$1] (col1 of file1) is being compared to a[$2] (col2 of file 2)
# {a[$1]=$2;next}: in file1, use col1 (headers) to remember col2 (ORF1_)
# if match found in file2,
# then print all of file2 plus append $2 from file1
# if match not found, print file2 as it is

# likewise, add annotation for the 'probable ORF1'
cat "$1"_probableORF1.txt \
| awk -F "_-_ORF" '{print $1 "\t" "probORF1_"}' \
| awk 'BEGIN{print "#comment"}{print}END{print ""}' \
| sed '/^$/d' \
> "$1"_probORF1_extendedheaders.tmp 

awk 'FNR==NR {a[$1]=$2;next} {if ($2 in a) print $0,a[$2]; else print $0;}' \
"$1"_probORF1_extendedheaders.tmp "$1"_annotated_with_ORF1.tmp \
> "$1"_annotated_with_ORF1_probORF1.tmp

# replace any cases of "ORF1_probORF1_" with just "ORF1"
cat "$1"_annotated_with_ORF1_probORF1.tmp \
| awk '{print $1 "\t" $2 "\t" $3$4}' \
| sed 's/ORF1_probORF1_/ORF1_/g' \
> "$1"_annotated_with_ORF1_probORF1_filt.tmp

# then add annotation for ORF2
cat "$1"_confirmedORF2.txt \
| awk -F "_-_ORF" '{print $1 "\t" "ORF2_"}' \
| awk 'BEGIN{print "#comment"}{print}END{print ""}' \
| sed '/^$/d' \
> "$1"_confirmedORF2_extendedheaders.tmp 

awk 'FNR==NR {a[$1]=$2;next} {if ($2 in a) print $0,a[$2]; else print $0;}' \
"$1"_confirmedORF2_extendedheaders.tmp "$1"_annotated_with_ORF1_probORF1_filt.tmp \
> "$1"_annotated_with_ORF1_probORF1_ORF2.tmp

# reorder the columns
cat "$1"_annotated_with_ORF1_probORF1_ORF2.tmp \
| awk '{print $3$4 "\t" $1 "\t" $2}' \
> "$1"_original_extended_with_orf_annotations.txt

# make a temp version 
# (only the original headers with ORF labels)
cat "$1"_annotated_with_ORF1_probORF1_ORF2.tmp \
| awk '{print $3$4 "\t" $1}' \
> "$1"_original_with_orf.tmp

# extract fasta and rename original headers with ORF annotation
# ORF label in col1; original header name in col2
awk 'NR==FNR{a[">"$2]=$1;next} {sub(/^>/,"&"a[$0])} 1' \
"$1"_original_with_orf.tmp /Volumes/Backup/Results/$2/Verified/"$1"_L1_verified.fasta \
> "$1"_L1_annotated.fasta

# move to a more convenient location
mv "$1"_L1_annotated.fasta /Volumes/Backup/Annotated_L1s/$2/

# remove temporary files
rm *.tmp 

# count the number of L1s that are :
# both orf (>ORF1_ORF2_* or >probORF1_ORF2_*),  
# orf1-only (>ORF1_* or >probORF1_*), and
# orf2-only (>ORF2_*)
cd /Volumes/Backup/Annotated_L1s/$2/

cat "$1"_L1_annotated.fasta \
| grep ">" \
| wc -l \
| awk '{print "total number of L1s: " $1}' \
> "$1"_L1_orf_counts.txt

grep -c "ORF1_ORF2_" "$1"_L1_annotated.fasta \
| awk '{print "# both orf L1s (>ORF1_ORF2_* or >probORF1_ORF2_*): " $1}' \
>> "$1"_L1_orf_counts.txt

cat "$1"_L1_annotated.fasta \
| grep "ORF1" \
| grep -v "ORF2" \
| wc -l \
| awk '{print "# orf1-only L1s (>ORF1_* or >probORF1_*): " $1}' \
>> "$1"_L1_orf_counts.txt

cat "$1"_L1_annotated.fasta \
| grep "ORF" \
| grep -v "ORF1" \
| wc -l \
| awk '{print "# orf2-only L1s (>ORF2_*): " $1}' \
>> "$1"_L1_orf_counts.txt

cat "$1"_L1_annotated.fasta \
| grep ">" \
| grep -v "ORF" \
| wc -l \
| awk '{print "# no orf L1s (>*): " $1}' \
