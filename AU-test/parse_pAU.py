#!/usr/bin/python


import re
from sys import argv


infile = argv[1]
#  'PairedSitesSupport\->{.*,\K[0-9]\.[0-9]*},Likel' 

with open("possible_trees") as fin:
	top_count = len(fin.readlines())

print "id\ttopology\tpAU"
with open(infile) as f:
	i = 1
	for line in f.readlines():
		m = re.search('PairedSitesSupport->{([0-9]\.[0-9]*,){5}([0-9]\.[0-9]*)}', line)
		p = re.search('Phylogeny->{(.*)},Subs', line)
		if m:
			print "top%i\t" %i,
			print p.groups()[0] + "\t",
			print(m.groups()[1] + "\n"),
			i += 1 
			if i % (top_count +1)  == 0:
				i = 1
		

