#!/bin/bash

indir=$1

# starting tree estiamtion in parallel
ls $indir/*.fa \
	| xargs -P8 -I {} ./gf-iqtree.sh {} &&

# collect gene trees
cat $(find $indir -name "*.treefile") > phylo-gftrees.tree
# collect AU values
echo -e "fragment_fname\ttopology\tpAU" > phylo-gftrees.au
cat $(find $indir -name "*.au") >> phylo-gftrees.au



