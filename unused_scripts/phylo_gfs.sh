#!/bin/bash

# Usage: ./phylo_gfs.sh <input_dir>

# Requirements: 'phylo_gfs.sh', 'genome_fragments.py', 'check_gfs.py',
#               'iqtree_parallel.sh', 'iqtree_au_test.sh', and 'parse_iqlog.py'
#               must be in the same directory.

# The following tools must be defined in $PATH:
# python 3
# parallel

# find input FASTAs
FASTAS=$(find $1 -name '*.clean.concat.fa' -printf '-f %p ')
MIN=25000  # minimum fragment size
MAX=250000  # maximum fragment size
STEP=25000  # step of size increase in range of fragment sizes
N_SAMPLES=500  # number of random fragments to sample

# generate random genome fragments (GF) of different sizes
for SIZE in $(seq $MIN $STEP $MAX); do
  # append command to joblist
  echo "python genome_fragments.py $FASTAS -c -n 0.2 -s $SIZE -r $N_SAMPLES" >> genome_fragments.jobs
done
# run multiple instances of 'genome_fragments.py' in parallel
# WARNING: memory intensive!
# ~ 27 GB RAM per instance for ~ 570 Mb alignment of 51 samples
# WARNING: very slow for large fragment sizes!
# ~ 27 days for generating 500 random 250 Kb fragments for an alignment of 51 samples
cat genome_fragments.jobs | parallel -j 5 &&
# remove intermediate files
rm genome_fragments.jobs

# rename FASTA headers
sed -Ei 's/\sScaffold_.*//g' *.fa

# create file with alternative tree topologies in newick format
cat <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),WOAK);")\
    <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),WOAK);")\
    <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),WOAK);")\
    <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),WOAK);")\
    <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),WOAK);")\
    <(echo "(((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),WOAK);")\
    <(echo "(((((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),WOAK);")\
    <(echo "(((((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),WOAK);")\
    <(echo "(((((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),WOAK);")\
    <(echo "(((((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),WOAK);")\
    <(echo "(((((MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),WOAK);")\
    <(echo "(((((MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),(WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24)),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),WOAK);")\
    <(echo "((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08)),((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110))),WOAK);")\
    <(echo "((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14)),((RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110))),WOAK);")\
    <(echo "((((WA720,WA733,WA746,WA806,WA808,GNP01,GNP04,GNP05,ZNP01,SNR2,ETH1,ETH2,ETH3,MF06,MF22,MF24),(MTNP09,BNP02,SUN3,KKR01,KKR08,V23,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110)),((LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,MA1,SGR01,SGR05,SGR07,SGR13,SGR14),(RET1,RET3,RET4,RET5,RET6,RETRot1,RETRot2,ISC01,ISC04,ISC08))),WOAK);")\
    > topologies.tree

for SIZE in $(seq $MIN $STEP $MAX); do
  # create GFs directory and move GF FASTAs to it
  mkdir GFs_${SIZE}bp && mv GF*_${SIZE}bp_*.fa GFs_${SIZE}bp
  # copy scripts to GF directory
  cp check_gfs.py iqtree_parallel.sh iqtree_au_test.sh parse_iqlog.py GFs_${SIZE}bp
  # change into GF directory
  cd GFs_${SIZE}bp
  # calculate proportion of N's per sequence per GF and save separate lists of good and bad GFs
  for GF in $(ls -v *.fa); do
    echo "sed '/^[^>]/ s/[^N]//gi; /^\s*$/d' $GF | python check_gfs.py $GF $SIZE 51 > ${GF/.fa/.percent_n}" >> check_gfs_${SIZE}bp.jobs
  done
  cat check_gfs_${SIZE}bp.jobs | parallel -j 20 && rm check_gfs_${SIZE}bp.jobs
  # randomly sample 100 GFs from the list of good GFs
  shuf -n 250 -o au_test.gfs good.gfs
  # change to parent directory
  cd ..
done

# perform AU test on each GF with IQ-TREE
for SIZE in $(seq $MIN $STEP $MAX); do
  echo "cd GFs_${SIZE}bp && bash ./iqtree_parallel.sh . 10 $SIZE" >> iqtree_parallel.jobs
done
# run multiple instances of 'iqtree_parallel.sh' in parallel
# WARNING: each job will use 10 CPUs!
cat iqtree_parallel.jobs | parallel -j 5 &&
# remove intermediate files
rm iqtree_parallel.jobs

# combine all 'phylo_GFs_*bp.au' files
cat $(ls -v *.au) | sed -r '2,$s/Fragment\tTopology\tpAU//g; /^\s*$/d' > combined.au

# run Astral-III with gene trees generated from best GF size
ASTRAL=/home/rcoimbra/software/astral/astral.5.6.3.jar
java -jar $ASTRAL -i phylo-gftrees.tree -o astral.tree

# create NEXUS input for PhyloNet
cat <(echo -e "#NEXUS\n\nBEGIN TREES;\n") \
    <(cat -n phylo-gftrees.tree) \
    <(echo -e "\nEND;\n\nBEGIN PHYLONET;\n\nMCMC_GT (gt1-gt100);\n\nEND;") \
  | sed -r "s/^\s+([0-9]+)\t/Tree gt\1 = /" > phylo-gftrees.nex

# run PhyloNet
java -jar /home/rcoimbra/software/PhyloNet_3.7.1.jar phylo-gftrees.nex > phylonet.out


# number of good GFs
for d in GFs_*; do cd $d; (for f in $(ls -v *.percent_n); do python test.py $f; done) | wc -l; cd ..; done
# number of scaffolds represented among GFs
for d in GFs_*; do cd $d; (for f in $(ls -v *.percent_n); do python test.py $f; done) | grep -Po 'Scaffold_\K\d+' | sort -n | uniq | wc -l; cd ..; done