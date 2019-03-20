SHELL := /bin/bash

MIN=1000
MAX=100000
STEP=1000
SAMPLE_SIZE=100

all: \
	$(subst %.pAU.tsv, %.pAU.out, $(wildcard *.pAU.out)) \
	size_dependent_pAU.final.txt  

#	*.pAU.out \

Bac00_scaffold*.fasta: steps
	python get_genomic_fragments_random.py \
		-i Bac00.sorted.mkdup.consensus.fasta \
		-i Bbo02.sorted.mkdup.consensus.fasta \
		-i Bph03.sorted.mkdup.consensus.fasta \
		-i Mno00.sorted.mkdup.consensus.fasta \
		-i Bbo01.sorted.mkdup.consensus.fasta \
		-i Bmu.merged.consensus.fasta \
		-i Egl00.sorted.mkdup.consensus.fasta \
		-i ../correct_repeatmasking/balMys.scaf.masked.filtered_100kb.fa \
		-s $(STEP) -r $(SAMPLE_SIZE) -c 1 --min $(MIN) --max $(MAX) 


*.pAU.tsv: Bac00_scaffold*.fasta
	for step in $$(seq $(MIN) $(STEP) $(MAX)); \
	do \
		ls Bac00_scaffold_*_$$step_*.fasta > fragment_filelist_$$step ;\
		./test_topologies.sh fragment_filelist_$$step fragment_filelist_$$step.pAU.tsv;\
	done



size_dependent_pAU.final.txt: *.pAU.tsv
	echo -e 'id\ttree\tfrag_size\tpAU' > $@
	cat $^ | grep -v "^id" >> $@

