#echo ">seq1" > seqA.fa; for i in $(seq 1 10000); do shuf -n 1 bases | tr -d "\n" >> seqA.fa; done;
#echo "" >> seqA.fa

function mod_sequence {
    cp $1 $2; # copy first sequence
    for i in $(seq 1 1); do 
        # creates random dinucleotide to search and replace in order to modify the sequence
        x="$(shuf -n 2 bases | tr -d '\n')"; # search sequence
        y="$(shuf -n 2 bases | tr -d '\n')"; # replace sequence
        sed -i "s/$x/$y/g" $2; # perform substitution
    done
}

# apply the randomozaion
mod_sequence seqA.fa seqB.fa;
mod_sequence seqA.fa seqC.fa;
mod_sequence seqA.fa seqD.fa;

# find ambiguos sites
for f in $(ls seq?.fa); 
    do     ./../find_sites_in_fasta.py -i $f -s "N" -s "Y" -s "R" -s "W" -s "K" -s "M" -s "S" -o ${f/.fa/.exclude.bed};
done

# merge bed files
cat *.exclude.bed  | bedtools sort -i - | bedtools merge -d 0 -c 4 -o collapse > all.exclude.bed
# create complement of bed files
# make sure to create a .genome file
bedtools complement -i all.exclude.bed -g seq1.genome > all.include.bed

# create genome fragments
python ../get_genomic_fragments_v2.py  -i seqA.fa \
                                    -i seqB.fa \
                                    -i seqC.fa \
                                    -i seqD.fa \
                                    -s 500 \
                                    -b all.include.bed\
                                    -c 
