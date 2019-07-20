'''

'''

from Bio import SeqIO
import argparse
from collections import defaultdict
from itertools import chain
from operator import itemgetter



#out_dic = defaultdict(list)


def find_char(id,s, ch):
    l = [i for i, ltr in enumerate(s) if ltr == ch]
    #print list(l)
    N_regions = []
    for g in grouper(l):
#      print minmax(g)
        N_regions.append(minmax(g))
    return N_regions 

def grouper(iterable):
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= 1:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group
        
def minmax(iterable):
    return [min(iterable),max(iterable)+1]



def main():
	seqs = SeqIO.index(options.input_fasta,"fasta")
	sites_to_find = [e.strip() for e in options.sites]
	with open(options.out_bed,"w") as fout:
		for record in seqs:
			record = seqs.get(record)
			seq = record.seq
			list_per_record = []
			for char in sites_to_find:
				for e in find_char(record.id,seq, char):
			    		list_per_record.append([record.id] + e +  [char])
			for l in sorted(list_per_record, key=itemgetter(2)):
			    fout.write("\t".join([str(e) for e in l])+"\n")			    
	return True
            
            
if __name__ == '__main__':
    program_description = """
    ### This is find_sites_in_FASTA.py ###

    (c) Fritjof Lammers 2016

	This script loads a FASTA and looks for specified base(e.g. ambiguities, Ns, etc) and outputs a BED file indicating the coordinates of 	identified sites. 

    """
    print program_description

    parser = argparse.ArgumentParser(description="Genomic fragments generator.")
    parser.add_argument('-i', '--input_fasta', required=True, help='Input file (FASTA)')
    parser.add_argument('-s', '--sites', action="append", required=True, help='Fragment size')
    parser.add_argument('-o', '--out_bed', required=True, help='Output file (BED)')

    options = parser.parse_args()


    main()
