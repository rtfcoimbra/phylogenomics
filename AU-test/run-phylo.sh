#!/bin/bash

indir=$1
threads=8
pattern=".fa"

# starting tree estiamtion in parallel
echo "Starting tree reconstruction"
echo "for $(ls $indir/*${pattern} | wc -l) genome fragments."
echo -e "doing ${threads} in parallel... (this might still take a while)...\n"
ls $indir/*${pattern} \
	| xargs -P${threads} -I {} ./gf-iqtree.sh {} &&

echo -e "\ndone!"

echo "Collecting data"
# collect gene trees
cat $(find $indir -name "*.treefile") > phylo-gftrees.tree

# collect AU values
echo -e "fragment_fname\ttopology\tpAU" > phylo-gftrees.au
cat $(find $indir -name "*.au") >> phylo-gftrees.au

echo "Done. Finished"

