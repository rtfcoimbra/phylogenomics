from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Nexus import Nexus
from os.path import basename
from sys import argv

fin = argv[1]
fout = basename(fin).split('.')[0] + '.nex'

alignment = AlignIO.read(fin, 'fasta')
minimal_record = '#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=dna; end;'
n = Nexus.Nexus(minimal_record)
n.alphabet = IUPAC.ambiguous_dna
for record in alignment:
    n.add_sequence(record.id, str(record.seq))
n.write_nexus_data(fout, interleave=False)
