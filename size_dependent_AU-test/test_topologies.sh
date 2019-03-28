#!/bin/bash

FILE_LIST=$1
OUT=$2

TEST_TREES_FILE=possible_trees

echo "" > $OUT
rm  ${FRAGMENT_FILE}.trees
mkdir "${FRAGMENT_FILE}_data"

while read FRAGMENT_FILE;
do

sed -i 's/_.*//g' $FRAGMENT_FILE
echo -e "TestTopologies[\"$FRAGMENT_FILE\",\"${TEST_TREES_FILE}\",SubstitutionModel->GTR[Optimum,Empirical]:GI[Optimum]:4,NReplicates->1000],\"${FRAGMENT_FILE/.fasta/.test_topo}\", SaveReport" > $INFILES$

## TREE GENERATION 
raxmlHPC-SSE3 -m GTRGAMMAI -n ${FRAGMENT_FILE}.raxml -s $FRAGMENT_FILE -p 432434 # create ML tree
echo -e "$FRAGMENT_FILE\t $(cat RAxML_bestTree.$FRAGMENT_FILE.raxml)" >>  trees
tree_estimate=$(sed 's/\:[0-9][\.]*[a-z0-9]*[-]*[0-9]*//g; s/\s//g' RAxML_bestTree.${FRAGMENT_FILE}.raxml)
rm RAxML_*.${FRAGMENT_FILE}*
if [ "${tree_estimate}" == "$(cat best_tree_raxml)" ];
then
	echo "[$0 / ${FRAGMENT_FILE} ]!!! Trees match $line !!!"
	# STATISTICAL LIKELIHOOD ESTIMATION
	raxmlHPC-SSE3 -m GTRGAMMAI -n ${FRAGMENT_FILE}.raxml -s ${FRAGMENT_FILE}.fragment_file.fa -f g -p 432434 -z $TEST_TREES_FILE # get Likelihood stats
	mv  RAxML_perSiteLLs.${FRAGMENT_FILE}.raxml  RAxML_perSiteLLs.${FRAGMENT_FILE}.sitelh
	makermt -b 0.01 --puzzle RAxML_perSiteLLs.${FRAGMENT_FILE}.sitelh && consel  RAxML_perSiteLLs.${FRAGMENT_FILE} && catpv  RAxML_perSiteLLs.${FRAGMENT_FILE} | \
	grep "[0-9].*" | awk '{print "top"$3, $5}' | sort -V -k 1,1 >> $OUT
	rm RAxML_*.${FRAGMENT_FILE}*
else
	echo "[$0 / ${FRAGMENT_FILE} ] !!! for ${line} tree did not match best tree !!!"
fi

done < $FILE_LIST

