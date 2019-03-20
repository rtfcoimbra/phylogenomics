#!/bin/bash

echo "... in progress: tree and AU calucation for $1"
# set path to alternative trees file, if none given AU test will be skipped
alt_trees=example/example.alt_trees
# run iqtree with (-au) option 10000 replicates and ultra-fast bootstrap (1000 replicates)
iqtree -s $1 -z ${alt_trees} -bb 1000  -zb 10000 -au > gf-iqtree.log &&
# extract AU values from iqtree logfile
python parse_iqlog.py $1.log > $1.au

