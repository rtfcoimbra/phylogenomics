SHELL := /bin/bash

GENOME_FILE=
EXCLUDE_FILES=$(patsubst %.cosensus.fasta,%.exclude_sites.bed, $(INPUT_FILES)

# Check that files exist
sites: \
	$(patsubst %.cosensus.fasta,%.exclude_sites.bed, $(INPUT_FILES)

coordinates:
	all_genomes.exclude_sites.bed	

# Identify ambiguos and gap sites
%.exclude_sites.bed: %.consensus.fasta
	./find_sites_in_fasta.py -i $< -s "N" -s "Y" -s "R" -s "W" -s "K" -s "M" -s "S" -o $@

# merge exclude sites to a common bed file
all_genomes.exclude_sites.bed: $(EXCLUDE_FILES)
	bedtools unionbedg -i $(EXCLUDE_FILES) | bedtools merge -d 0 -c 4 -o distinct > $@; \

# create include file from exclude file (complement)
all_genomes.include_sites.bed: all_genomes.exclude_sites.bed
	bedtools complement -i all_genomes.exclude_sites.bed -g $(GENOME_FILE)

# create genome fragments
gfs_done: all_genomes.include_sites.bed
	python get_genomic_fragments_v2b.py \
			-i scaffold1.fasta \
			-i scaffold1_b.fasta \
			-s $(FRAGMENT_SIZE)
			-b $< \
			-c


