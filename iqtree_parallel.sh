#!/bin/bash

# Usage: ./iqtree_parallel.sh <input_dir> <n_jobs> <fragment_size>

# Requirements: 'phylo_GFs.sh', 'genome_fragments.py', 'iqtree_parallel.sh',
#               'iqtree_AU_test.sh', and 'parse_iqlog.py' must be in the same
#               directory.

# The following tools must be defined in $PATH:
# python 3
# parallel
# iqtree

echo "Starting phylogenetic analyses for $(ls $1/*.fa | wc -l) genome fragments."
echo -e "Running $2 instances in parallel.\n"

# run multiple instances of IQ-TREE and parse_iqlog.py in parallel
ls -v $1/*.fa | parallel -k -j $2 'bash ./iqtree_au_test.sh {}' &&

echo -e "Done.\n"
echo "Now collecting data..."

# collect gene trees
cat $(find $1 -name "*.treefile") > phylo_GFs_${3}bp.tree
# collect p-AU values
echo -e "Fragment\tTopology\tpAU" > phylo_GFs_${3}bp.au
cat $(find $1 -name "GF*.au") >> phylo_GFs_${3}bp.au
# edit fragment name in the first column of the '.au' file
sed -ri "s/^.*_($3)bp_.*log/\1/" phylo_GFs_${3}bp.au

# copy outputs to parent directory
cp phylo_GFs_${3}bp.tree phylo_GFs_${3}bp.au ..

echo "Done."
