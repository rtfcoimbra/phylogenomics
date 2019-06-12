#!/bin/bash

# Usage: ./clean_masked_fasta.sh <input_dir> <output_dir> <n_cpus>

# Requirements: 'clean_masked_fasta.sh' and 'concatenate_fasta.py' must be in the same directory.

# The following tools must be defined in $PATH:
# python 3
# parallel
# bedtools

# path to BED file of non-repetitive regions of scaffolds >= 1 Mb
BED=/gendata_aj/rcoimbra/reference_assembly/kordofan_giraffe/no_repeats_1mb_scaffolds.bed

# extract non-repetitive regions from masked FASTA
for FASTA in $1/*.fa; do
  CLEAN_FASTA=$2/$(basename ${FASTA/.fa/.clean.fa})
  echo "bedtools getfasta -fi $FASTA -bed $BED -fo $CLEAN_FASTA" >> $2/getfasta.jobs
done
cat $2/getfasta.jobs | parallel -j $3

# concatenate FASTA regions per sample
for CLEAN_FASTA in $2/*.clean.fa; do
  CONCAT_FASTA=$2/$(basename ${CLEAN_FASTA/.fa/.concat.fa})
  echo "python ./concatenate_fasta.py $CLEAN_FASTA $CONCAT_FASTA" >> $2/concatenate_fasta.jobs
done
cat $2/concatenate_fasta.jobs | parallel -j $3

# mask remaining ambiguous sites
#sed -i '/^[^>]/ s/[^AGTC]/N/gi' $2/*.clean.concat.fa

# index FASTAs
for CONCAT_FASTA in $2/*.clean.concat.fa; do
  echo "samtools faidx $CONCAT_FASTA" >> $2/samtools_faidx.jobs
done
cat $2/samtools_faidx.jobs | parallel -j $3

# calculate proportion of N's per FASTA
for CONCAT_FASTA in $1/*.clean.concat.fa; do
  BASE_COUNT=$(grep -v '^>' $CONCAT_FASTA | tr -d '\n' | wc -c)
  N_COUNT=$(grep -v '^>' $CONCAT_FASTA | tr -cd N | wc -c)
  N_PERCENT=$(python -c "print(f'{round(($N_COUNT / $BASE_COUNT) * 100, 2)}')")
  echo -e "$(basename $CONCAT_FASTA)\t$BASE_COUNT\t$N_COUNT\t$N_PERCENT" >> $2/n_percent.tmp
done
cat <(echo -e "FILE\tBASE_COUNT\tN_COUNT\tN_PERCENT") <(sort -gr -k4,4 $2/n_percent.tmp) > $2/n_percent.tbl

# remove intermediate files
rm $2/*.jobs $2/*.clean.fa $2/n_percent.tmp
