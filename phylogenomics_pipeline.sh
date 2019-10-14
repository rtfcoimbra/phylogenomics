# The following scripts must be in the same directory:
# 'giraffe_au_test.sh'
# 'genome_fragments.py'
# 'check_gfs.py'
# 'iqtree_parallel.sh'
# 'iqtree_au_test.sh'
# 'parse_iqlog.py'

# The following tools must be defined in $PATH:
# python 3
# parallel


################################################################################
#                      GF size-dependent AU topology test                      #
################################################################################

# randomly sample two representatives per giraffe (sub)species
# only quality filtered genome FASTAs with <= 20% of N's were considered
shuf -e -n 2 WA720 WA733 WA746 WA806 WA808
shuf -e -n 2 GNP01 GNP04 GNP05 SNR2 ZNP01
shuf -e -n 2 ETH1 ETH2 ETH3 MF06 MF22 MF24
shuf -e -n 2 ISC04 ISC08 RET1 RET3 RET4 RET5 RET6 RETRot1 RETRot2
shuf -e -n 2 MA1 SGR01 SGR05 SGR07 SGR13 SGR14
shuf -e -n 2 BNP02 KKR01 KKR08 MTNP09 SUN3 V23
shuf -e -n 2 ENP11 ENP16 ENP19 ENP20 HNB102 HNB110

# directory containing input FASTAs
DIR=/home/rcoimbra/results/angsd_dofasta
# FASTAs of randomly sampled giraffe (sub)species representatives
FASTAS="-f $DIR/WA733.clean.concat.fa -f $DIR/WA808.clean.concat.fa -f $DIR/GNP01.clean.concat.fa -f $DIR/SNR2.clean.concat.fa -f $DIR/ETH2.clean.concat.fa -f $DIR/MF06.clean.concat.fa -f $DIR/ISC04.clean.concat.fa -f $DIR/RETRot2.clean.concat.fa -f $DIR/SGR05.clean.concat.fa -f $DIR/SGR14.clean.concat.fa -f $DIR/KKR08.clean.concat.fa -f $DIR/V23.clean.concat.fa -f $DIR/ENP16.clean.concat.fa -f $DIR/ENP19.clean.concat.fa -f $DIR/WOAK.clean.concat.fa"
MIN=50000  # minimum fragment size
MAX=600000  # maximum fragment size
STEP=50000  # step of size increase in range of fragment sizes
N_SAMPLES=200  # number of random fragments to sample

# generate random GFs of different sizes
for SIZE in $(seq $MIN $STEP $MAX); do
  # append command to a joblist
  echo "python genome_fragments.py $FASTAS -c -n 0.2 -s $SIZE -r $N_SAMPLES" >> genome_fragments.jobs
done
# run multiple instances of 'genome_fragments.py' in parallel
# WARNING: memory intensive and very slow for large fragment sizes
cat genome_fragments.jobs | parallel -j 10

# rename FASTA headers
sed -Ei 's/\sScaffold_.*//g' *.fa

# set up a variable for each species containing its representatives
NOR="WA733,WA808,GNP01,SNR2,ETH2,MF06"  # northern giraffe
RET="ISC04,RETRot2"  # reticulated giraffe
MAS="SGR05,SGR14"  # masai giraffe
SOU="KKR08,V23,ENP16,ENP19"  # southern giraffe

# create file with alternative tree topologies in newick format
cat <(echo "((((($NOR),($RET)),($MAS)),($SOU)),WOAK);")\
    <(echo "((((($NOR),($RET)),($SOU)),($MAS)),WOAK);")\
    <(echo "((((($NOR),($MAS)),($RET)),($SOU)),WOAK);")\
    <(echo "((((($NOR),($MAS)),($SOU)),($RET)),WOAK);")\
    <(echo "((((($NOR),($SOU)),($MAS)),($RET)),WOAK);")\
    <(echo "((((($NOR),($SOU)),($RET)),($MAS)),WOAK);")\
    <(echo "((((($MAS),($SOU)),($NOR)),($RET)),WOAK);")\
    <(echo "((((($MAS),($SOU)),($RET)),($NOR)),WOAK);")\
    <(echo "((((($MAS),($RET)),($NOR)),($SOU)),WOAK);")\
    <(echo "((((($MAS),($RET)),($SOU)),($NOR)),WOAK);")\
    <(echo "((((($SOU),($RET)),($MAS)),($NOR)),WOAK);")\
    <(echo "((((($SOU),($RET)),($NOR)),($MAS)),WOAK);")\
    <(echo "(((($NOR),($RET)),(($MAS),($SOU))),WOAK);")\
    <(echo "(((($NOR),($MAS)),(($RET),($SOU))),WOAK);")\
    <(echo "(((($NOR),($SOU)),(($MAS),($RET))),WOAK);")\
    > topologies.tree

# for each GF size
for SIZE in $(seq $MIN $STEP $MAX); do
  # create GFs directory and move GF FASTAs to it
  mkdir GFs_${SIZE}bp && mv GF*_${SIZE}bp_*.fa GFs_${SIZE}bp
  # copy scripts to GF directory
  cp check_gfs.py iqtree_parallel.sh iqtree_au_test.sh parse_iqlog.py GFs_${SIZE}bp
  # change into GF directory
  cd GFs_${SIZE}bp
  # calculate proportion of N's per sequence per GF and save separate lists of good and bad GFs
  for GF in $(ls -v *.fa); do
    # append command to a joblist
    echo "sed '/^[^>]/ s/[^N]//gi; /^\s*$/d' $GF | python check_gfs.py $GF $SIZE 15 > ${GF/.fa/.percent_n}" >> check_gfs_${SIZE}bp.jobs
  done
  # run multiple instances of 'check_gfs.py' in parallel
  cat check_gfs_${SIZE}bp.jobs | parallel -j 10
  # change to parent directory
  cd ..
done

# perform AU test on each GF with IQ-TREE
for SIZE in $(seq $MIN $STEP $MAX); do
  # append command to a joblist
  echo "cd GFs_${SIZE}bp && bash ./iqtree_parallel.sh . 10 $SIZE" >> iqtree_parallel.jobs
done
# run multiple instances of 'iqtree_parallel.sh' in parallel
# WARNING: each job will use 10 CPUs
cat iqtree_parallel.jobs | parallel -j 5

# combine all 'phylo_GFs_*bp.au' files
cat $(ls -v *.au) | sed -r '2,$s/Fragment\tTopology\tpAU//g; /^\s*$/d' > combined.au


################################################################################
#                          Whole-genome phylogenomics                          #
################################################################################

# BED file of scaffolds >= 1 Mb
BED=/gendata_aj/rcoimbra/reference_assembly/kordofan_giraffe/kordofan-giraffe_1mb_scaffolds.bed
# sort BED by scaffold size and split it into new BEDs with up to 10 scaffolds
cat $BED | sort -Vk 3 | split -l 1 -a 3 -d - kordofan-giraffe_1mb_scaffolds.bed.
# directory containing input FASTAs
DIR=/home/rcoimbra/results/angsd_dofasta
# find input FASTAs (only quality filtered genome FASTAs with <= 20% of N's were considered)
FASTAS=$(find $DIR -name '*.clean.concat.fa' ! -name 'ISC01.clean.concat.fa' ! -name 'LVNP*.clean.concat.fa' -printf '-f %p ')
# generate non-overlapping GFs of 450 Kb
ls -v kordofan-giraffe_1mb_scaffolds.bed.* | parallel -j 10 "python genome_fragments.py $FASTAS -b {} -c -n 0.2 -s 450000"
# rename FASTA headers
sed -Ei 's/\sScaffold_.*//g' *.fa

Scaffold_6
Scaffold_60
Scaffold_4
Scaffold_2
Scaffold_1
Scaffold_3
Scaffold_0

# calculate proportion of N's per sequence per GF and save separate lists of good and bad GFs
for GF in $(ls -v *.fa); do
  # append command to a joblist
  echo "sed '/^[^>]/ s/[^N]//gi; /^\s*$/d' $GF | python check_gfs.py $GF 450000 $(grep -c '>' $GF) > ${GF/.fa/.percent_n}" >> check_gfs_450kb.jobs
done
# run multiple instances of 'check_gfs.py' in parallel
cat check_gfs_450kb.jobs | parallel -j 20

# perform ultrafast model selection followed by tree reconstruction with 1000 ultrafast bootstrap replicates
ls -v *.fa | parallel -j 30 "iqtree -s {} -o 'WOAK' -bb 1000"
# collect gene trees
cat $(find . -name '*.treefile') > GFs_450kbp.tree

# run Astral-III with gene trees generated from best GF size
ASTRAL=/home/rcoimbra/software/astral/astral.5.6.3.jar
java -jar $ASTRAL -i phylo-gftrees.tree -o astral.tree
