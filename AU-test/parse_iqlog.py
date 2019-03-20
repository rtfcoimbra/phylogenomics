#!/usr/env python

from sys import argv

infile = argv[1]
isAU = False
with open(infile) as fin:
    for line in fin.readlines():
        if line.startswith("TreeID") and not isAU:
            isAU = True
        elif isAU and line.startswith("Time"):
            isAU = False
            pass
        elif isAU:
            lsplit = line.split("\t")
            top = "top%s" % lsplit[0]
            pAU = lsplit[1]
            print("{infile}\t{top}\t{pAU}".format(infile=infile, top=top, pAU=pAU))


            # Output info
            # AU p-value
            # RSS residual sum of squares
            # d distance
            # c curvature
