#!/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sys import argv
import re

infile=argv[1]
outfile=argv[2]
seqs = SeqIO.index(infile, "fasta")
outseqs = []
chr_ids = [record.split(":")[0] for record in seqs.keys()]

for chr in sorted(set(chr_ids)):
    matching_chr = [seqid for seqid in seqs.keys() if seqid.startswith("%s:"%chr)]
    collect_sequences = [str(seqs[seqid].seq) for seqid in matching_chr]
    chrseq = "".join(collect_sequences)
    print(">%s" %chr)
#    print(chrseq)
    chrrecord = SeqRecord(Seq(chrseq, generic_dna),
                          id="%s" %chr,
                          description = "")
    outseqs.append(chrrecord)

SeqIO.write(outseqs, outfile, "fasta")




