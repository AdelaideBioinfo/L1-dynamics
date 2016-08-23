#!/bin/bash

# Re-screen unconfirmed ORF1 seqs against the library of probable L1 ORF1 domains

# $1 = species name
# $2 = order
# Example usage: probableORF1.sh Homo.sapiens Mammalia

# Count the 'known' domains from confirmed ORF1 seqs across species
#cat */*/*_ORF1_domains_known.txt \
#| awk '{print $2}' \
#| tr "," "\n" \
#| grep -v "^\s*$" \
#| sort \
#| uniq -c \
#| sort -bnr \
#> Known_domain_count.txt

# Count the 'unknown' domains from confirmed ORF1 seqs across species
#cat */*/*_ORF1_domains_unknown.txt \
#| awk '{print $2}' \
#| tr "," "\n" \
#| grep -v "^\s*$" \
#| sort \
#| uniq -c \
#| sort -bnr \
#> Unknown_domain_count.txt

# Find common domains between confirmed and unconfirmed ORF1 seqs
#cat Known_domain_count.txt | awk '{print $2}' > known_domains.txt
#cat Unknown_domain_count.txt | awk '{print $2}' > unknown_domains.txt
#comm -1 -2 <(sort known_domains.txt) <(sort unknown_domains.txt) > domain_library.txt
#3648 domains total

# Use these common domains to re-screen the unconfirmed ORF1 seqs into probable ORF1
# This is the library of probable L1 domains
# Compare this to the domains found in the unconfirmed seqs (recursively)

cd /media/Backup/ORF_validation/$2/$1

awk 'NR==FNR { a[$1] = $1; next} { for (k in a) if ($2 ~ a[k]) { print $0; break } }' \
../../domain_library.txt "$1"_ORF1_domains_unknown.txt > "$1"_ORF1_domains_probable.txt

# Extract protein seqs of the probable ORF1 seqs
cat "$1"_ORF1_domains_probable.txt \
| awk '{print $1}' \
| sort -u \
> "$1"_probableORF1.txt

cat "$1"_extended_ORF1_candidates_translation.fasta \
| awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' \
| awk -F"\t" 'BEGIN{while((getline k < "'$1'_probableORF1.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' \
> "$1"_probableORF1.fasta
