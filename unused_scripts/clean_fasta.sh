#!/bin/bash

# Usage: ./clean_fasta.sh <input_dir> <output_dir> <n_cpus>

# Requirements: `clean_fasta.sh`, `find_sites_in_fasta.py`, and `concat_fasta.py` must be in the same directory.

# The following tools must be defined in $PATH:
# python 3
# parallel
# bedtools

# create a BED file with the positions of ambiguous sites and N's
for FASTA in $1/*.fa; do
  EX_BED=$2/$(basename ${FASTA/.fa/.exclude.bed})
  echo "python ./find_sites_in_fasta.py -i $FASTA -s 'N' -s 'Y' -s 'R' -s 'W' -s 'K' -s 'M' -s 'S' -o $EX_BED" >> $2/find_sites.jobs
done
cat $2/find_sites.jobs | parallel -j $3

# merge all BED files
cat $2/*.exclude.bed | sort -V | bedtools merge -d 0 -c 4 -o collapse > $2/all.exclude.bed

# generate a genome `.tsv` file
cut -f 1,2 /gendata_aj/rcoimbra/reference_assembly/kordofan_giraffe/kordofan-giraffe_1mb_scaffolds.masked.fa.fai > $2/scaffolds.1mb.tsv

# create a BED file containing regions with no ambiguous sites or N's
bedtools complement -i $2/all.exclude.bed -g $2/scaffolds.1mb.tsv > $2/all.include.bed

# extract FASTA regions with no ambiguous sites or N's
for FASTA in $1/*.fa; do
  CLEAN_FASTA=$2/$(basename ${FASTA/.fa/.clean.fa})
  echo "bedtools getfasta -fi $FASTA -bed $2/all.include.bed -fo $CLEAN_FASTA" >> $2/getfasta.jobs
done
cat $2/getfasta.jobs | parallel -j $3

# concatenate FASTA regions per sample
for CLEAN_FASTA in $2/*.clean.fa; do
  CONCAT_FASTA=$2/$(basename ${CLEAN_FASTA/.fa/.concat.fa})
  echo "python ./concat_fasta.py $CLEAN_FASTA $CONCAT_FASTA" >> $2/concat_fasta.jobs
done
cat $2/concat_fasta.jobs | parallel -j $3

# remove intermediate files
rm $2/*.jobs
