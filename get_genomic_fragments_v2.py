#!/usr/env python

__author__ = 'Fritjof Lammers'
"""

Cut Genome alignment in fragments after filtering out masked sequence.

Time estimate:
8 x 600kb = 4.6 Mb
Fragment Size: 10kb
2 min @ local machine

whole genome:
8 x 2.3 Gb = 20.4 Gb



"""
import argparse
from Bio import AlignIO
from Bio import Align
from Bio import SeqIO
# from Bio import Seq
import sys
# from pybedtools import BedTool
from os.path import basename, dirname, isfile
import itertools


def load_BED(fname):
    '''
    This functions loads a BED file to generate a query list from chromosome names.
    Note that positions (columns 2 and 3) are not used!

    :param  fname: filename of BED file with chromosome sizes
    :return: list of chromosomes from BED file
    '''
    bed_list = []
    with open(fname) as fin:
        for line in fin.readlines():
            entry = [e.strip() for e in line.split("\t")]
            bed_list.append([entry[0], int(entry[1]), int(entry[2])])

    return bed_list



def clean_scaffold(alignment, sites_to_include):
    '''
    This function iterates over the alignment column by column.
    Builds up new alignment containing only pure ATCG columns.
    :param alignment:
    :return:
    '''
    print "clean scaffold %s" % alignment[0].id
    intervals_to_include = [e[1:3] for e in sites_to_include]
    print max([e[1] for e in intervals_to_include])
    print len(alignment[1])

    for i in itertools.chain.from_iterable([range(x,y) for x,y in intervals_to_include]):
        if i >= len(alignment[1]):
            break # don't travel longer than alignment
        clean_sites = alignment[:, i:i + 1]
        if "cleaned_alignment" in locals():
            cleaned_alignment += clean_sites  # add column
            if len(cleaned_alignment[0]) % 10000 == 0:
               sys.stdout.write(".")
               sys.stdout.flush()
            else:
               pass
        else:
            cleaned_alignment = clean_sites
            sys.stdout.write("add sites\n")

    try:
        return cleaned_alignment
    except:
        return False


def chop_scaffold(alignment, frag_size):
    '''
    This function chops the alignment into fragments of frag_size.
    If a fragment becomes smaller than frag_size it is discarded.
    :param alignment: Bio.Align object
    :param frag_size: Size of fragments in bp, (int)
    :return: List of fragment_name, fragment tuples
    '''

    fragment_list = []
    fragment_no = 0
    for i in xrange(0, len(alignment[0]), frag_size):
        fragment = alignment[:, i:i + frag_size]

        if len(fragment[0]) < frag_size:
            print "Fragment smaller then requested, skipping (%i < %i)" % (len(fragment[0]), frag_size)
            continue
        fragment_list.append(("%s_%04d" % (alignment[0].id, fragment_no), fragment))
        fragment_no += 1

    return fragment_list


def main():

    global f_ab, f_extra, bt_positions
    seqs = {}
    records = []
    fname_list = [basename(fpath) for fpath in options.input_files]
    include_sites = load_BED(options.bed_file) if options.bed_file else []
    exclude_list = load_BED(options.exclude_bed) if options.exclude_bed else []

    assert type(options.fragment_size == "int")

    for fpath in options.input_files:  # read input fasta files
        print "loading fasta files..."
        fname = basename(fpath)
        seqs[fname] = SeqIO.index(fpath, "fasta")  # load sequences of each fasta file as SeqIO.indexd
        records_per_fasta = seqs.get(fname).keys()  # create list of fasta headers (i.e. scaffolds)
        records.extend(records_per_fasta)
        print "..."+fname

    print "\n"

    records = set(records)  # make unique list of scaffolds, i.e. remove duplicates

    #bt_positions = BedToolPositions()

    for record in sorted(records):  # iterate over scaffolds
        sequences = []

        if options.exclude_bed and record in exclude_list:
            print "%s in exclude list, skipping" % record
            continue

        for fname in fname_list:  # iterate over fasta files loaded ...
            # print seq_key
            seq_to_add = seqs.get(fname).get(record)  # load scaffold from fasta file
            seq_to_add.id = fname.split(".")[0] + "_" + seq_to_add.id  # generate new header from filename
            seq_to_add.name = seq_to_add.id  # copy ID to name attribute
            sequences.append(seq_to_add)  # add sequence from the current scaffold

        min_alignment_length = min([len(sequence) for sequence in sequences])  # determine minimum alignment length

        per_chr_alignment = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in sequences]) \
            # create alignment and reduce scaffold length to shortest scaffold
        if options.clean:
            processed_alignment = clean_scaffold(per_chr_alignment, [entry for entry in include_sites if entry[0] == record]) # perform clean of aligned scaffold
            if not processed_alignment:
                print "Something went wrong while cleaning up alignment of scaffold %s" %record
                continue
        else:
            processed_alignment = per_chr_alignment

        fragments = chop_scaffold(processed_alignment, options.fragment_size) # create fragments
        for frag_id, frag_alignment in fragments: # for each fragment, save separate fasta file
            frag_fname = frag_id + ".fasta"
            with open(frag_id + ".fasta", "w") as fout:
                print "writing %s..." % frag_id,
                AlignIO.write(frag_alignment, fout, "fasta")
                print "done."
    return 1


if __name__ == '__main__':
    program_description = """
    ### This is Genomic Fragment Generator. ###

    (c) Fritjof Lammers 2015

    The program loads several (genomic) consenus sequences in FASTA format, aligns them by identical headers and saves
    alignments for fragments of given size. \n
    - Give each input FASTA with separate -i parameter, e.g. -i species1.fasta -i species2.fasta ... -i speciesN.fasta
    - Output files are written to working directory.
    - To restrict the program to specific regions (right now only chromosome/scaffold names considered), you can specifiy
    a BED file (optional).

    """
    print program_description

    parser = argparse.ArgumentParser(description="Genomic fragments generator.")
    parser.add_argument('-i', '--input_files', required=True, action="append", help='Input files')
    parser.add_argument('-s', '--fragment_size', type=int, default=100000, required=True, help='Fragment size')
    parser.add_argument('-b', '--bed_file', required=False, help='BED files with sites to include')
    parser.add_argument('-e', '--exclude_bed', required=False, help='Exclude bed file')
    parser.add_argument('-c', '--clean', required=False, help='Set if scaffolds should be cleaned', action="store_true")

    options = parser.parse_args()


    main()
