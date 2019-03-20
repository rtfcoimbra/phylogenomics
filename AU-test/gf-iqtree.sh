#!/bin/bash

# set path to alternative trees file
alt_trees=example/example.alt_trees

# run iqtree with (-au) option 10000 replicates
iqtree -s $1 -z ${alt_trees} -bb 1000  -zb 10000 -au -redo &&
# extract AU values from iqtree logfile
python parse_iqlog.py $1.log > $1.au
