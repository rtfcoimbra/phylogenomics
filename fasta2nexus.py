import os.path
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC

infile = sys.argv[1]
outfile = os.path.splitext(infile)[0] + '.nex'

with open(infile, 'rU') as fin:
    with open(outfile, 'w') as fout:
        count = SeqIO.convert(fin, 'fasta', fout, 'nexus', alphabet=IUPAC.ambiguous_dna)

print(f"Converted {count} records")
