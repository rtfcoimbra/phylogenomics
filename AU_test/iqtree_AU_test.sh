#!/bin/bash

# Usage: ./iqtree_AU_test.sh <input_file>

# Requirements: 'phylo_GFs.sh', 'genome_fragments.py', 'iqtree_parallel.sh',
#               'iqtree_AU_test.sh', and 'parse_iqlog.py' must be in the same
#               directory.

# The following tools must be defined in $PATH:
# python 3
# parallel
# iqtree

# path to file of tree topologies, if none is given AU test will be skipped
ALT_TREES=$(find .. -name 'topologies.tree' -printf '%p')

echo "Tree reconstruction for $(basename $1) in progress..."
echo "Approximately unbiased tree topology test for $(basename $1) in progress..."

# run IQ-TREE with ultrafast bootstrap, ultrafast model selection, and AU tree topology test
iqtree -s $1 -o 'WOAK' -bb 1000 -nm 10000 -z $ALT_TREES -zb 10000 -au &&

# extract p-AU values from IQ-TREE's log file
python parse_iqlog.py ${1}.log > ${1}.au
