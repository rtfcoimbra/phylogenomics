#!/bin/env python3

from sys import argv

infile = argv[1]
outfile = argv[2]
last_chr = ""

print("Starting to concating per chromosome / scaffold...")

with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin.readlines():
        if line.startswith(">"):
            chrom = line.strip().split(":")[0]
            if not last_chrom:
                fout.write("%s\n" % chrom)  # write header
                print(chrom.replace(">", "..."))  # write status to stdout
            elif chrom != last_chrom:
                fout.write("\n%s\n" % chrom)  # write header
                print(chrom.replace(">", "..."))  # write status to stdout
            else:
                pass
        else:
            fout.write(line.strip()) # write sequence
        last_chrom = chrom

print("Done.")
# EOF
