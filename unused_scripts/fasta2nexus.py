from Bio import SeqIO
from Bio.Alphabet import IUPAC
from os.path import basename
from sys import argv

infile = argv[1]
outfile = basename(infile).split('.')[0] + '.nex'

with open(infile, 'rU') as fin:
    with open(outfile, 'w') as fout:
        SeqIO.convert(fin, 'fasta', fout, 'nexus', alphabet=IUPAC.ambiguous_dna)
