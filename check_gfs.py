import sys
from Bio import SeqIO

filename = sys.argv[1]
frag_size = sys.argv[2]
n_seqs = sys.argv[3]
acceptable_seqs = 0
bad_seqs = 0

for seq_record in SeqIO.parse(sys.stdin, 'fasta'):
    sequence = str(seq_record.seq).upper()
    percent_n = round((sequence.count('N') / int(frag_size)) * 100, 2)
    print(f"{seq_record.id}\t{percent_n}%")
    if 20 < percent_n <= 50:
        acceptable_seqs += 1
    elif percent_n > 50:
        bad_seqs += 1

with open('good.gfs', 'a') as fout1, open('bad.gfs', 'a') as fout2:
    if bad_seqs == 0 and acceptable_seqs <= int(n_seqs) * 0.2:
        fout1.write(f"{filename}\t{acceptable_seqs}\t{bad_seqs}\n")
    else:
        fout2.write(f"{filename}\t{acceptable_seqs}\t{bad_seqs}\n")

    fout1.close()
    fout2.close()
