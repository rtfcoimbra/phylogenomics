#!/bin/bash

# Usage: ./iqtree_au_test.sh <input_file>

# Requirements: 'genome_fragments.py', 'iqtree_parallel.sh',
#               'iqtree_au_test.sh', and 'parse_iqlog.py' must be in the same
#               directory.

# The following tools must be defined in $PATH:
# python 3
# parallel
# iqtree

# path to file of tree topologies, if none is given topology test will be skipped
ALT_TREES=$(find .. -name 'topologies.tree' -printf '%p')

echo "Tree reconstruction for $(basename $1) in progress..."
echo "Approximately unbiased (AU) tree topology test for $(basename $1) in progress..."

# perform AU tree topology test with 10000 replicates using ultrafast model selection
iqtree -s $1 -o 'WOAK' -n 0 -z $ALT_TREES -zb 10000 -au &&

# extract p-AU values from IQ-TREE's log file
python parse_iqlog.py ${1}.log > ${1}.au
