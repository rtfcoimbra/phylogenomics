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
from os.path import basename, dirname
from os.path import basename, dirname
import itertools
import random

def contains_char(seq, aset):
    '''
    :param seq: string of sequence
    :param aset: set of character to be checked
    :return: True if character found, False if not
    '''
    for item in itertools.ifilter(aset.__contains__, seq):
        return True
    return False


def load_BED(fname):
    '''
    This functions loads a BED file to generate a query list from chromosome names.
    Note that positions (columns 2 and 3) are not used!

    :param  fname: filename of BED file with chromosome sizes
    :return: list of chromosomes from BED file
    '''
    chromosomes = []
    with open(fname) as fin:
        for line in fin.readlines():
            chromosomes.append(line.split("\t")[0])

    return chromosomes


def clean_scaffold(alignment):
    '''
    This function iterates over the alignment column by column.
    Builds up new alignment containing only pure ATCG columns.
    :param alignment:
    :return:
    '''
    print "clean scaffold %s" % alignment[0].id
    for i in xrange(0, len(alignment[1])):  # iterate over alignment
        sites = alignment[:, i:i + 1]
        badchar = False

        for nuc in sites:
            if str(nuc.seq).upper() in "NYRKMWSBDHV-":
                badchar = True
                break

        if not badchar:
            if "cleaned_alignment" in locals():
                cleaned_alignment += sites  # add column
                #if len(cleaned_alignment[0]) % 100 == 0:
                #    sys.stdout.write(".")
                #    sys.stdout.flush()
                #else:
                #    pass
            else:
                cleaned_alignment = sites
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
        fragment_list.append(("%s_%04d" % (alignment[0].id, fragment_no, frag_size), fragment))
        fragment_no += 1

    return fragment_list


def iterate_all(records,  fname_list, query_list, seqs):

    for record in records:  # iterate over scaffolds
      #        sequences = []
        if options.bed_file and record not in query_list :
            print "%s not in list, skipping" % record
            continue

        for fname in fname_list:  # iterate over fasta files loaded ...
            # print seq_key
            seq_to_add = seqs.get(fname).get(record)  # load scaffold from fasta file
            seq_to_add.id = fname.split(".")[0] + "_" + seq_to_add.id  # generate new header from filename
            seq_to_add.name = seq_to_add.id  # copy ID to name attribute
            sequences.append(seq_to_add)  # add sequence from the current scaffold

        min_alignment_length = min([len(sequence) for sequence in sequences])  # determine minimum alignment length

        per_chr_alignment = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in sequences])   # create alignment and reduce scaffold length to shortest scaffold
        if options.clean:
            processed_alignment = clean_scaffold(per_chr_alignment) # perform clean of aligned scaffold
            if not processed_alignment:
                print "Something went wrong while cleaning up alignment of scaffold %s" %record
                continue
        else:
            processed_alignment = per_chr_alignment

        fragments =  processed_alignment #fname_list, query_listchop_scaffold(processed_alignment, options.fragment_size) # create fragments

        for frag_id, frag_alignment in fragments: # for each fragment, save separate fasta file
            with open(frag_id + ".fasta", "w") as fout:
                print "writing %s..." % frag_id,
                AlignIO.write(frag_alignment, fout, "fasta")
                print "done."

        return fragments()

def iterate_random(records, fname_list, query_list, seqs, fragment_size):
    fragment_count = 0

    while fragment_count < options.random:
        record = random.choice(list(records))
        sequences = []
        if options.bed_file and record not in query_list :
            print "%s not in list, skipping" % record
            continue

        for fname in fname_list:  # iterate over fasta files loaded ...
            # print seq_key
            seq_to_add = seqs.get(fname).get(record)  # load scaffold from fasta file
            seq_to_add.id = fname.split(".")[0] + "_" + seq_to_add.id  # generate new header from filename
            seq_to_add.name = seq_to_add.id  # copy ID to name attribute
            sequences.append(seq_to_add)  # add sequence from the current scaffold

        min_alignment_length = min([len(sequence) for sequence in sequences])  # determine minimum alignment length

        per_chr_alignment = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in sequences])
        
        if min_alignment_length > (fragment_size*5):
            random_start = random.randint(0, min_alignment_length -(fragment_size*5))
            per_chr_alignment = per_chr_alignment[:, random_start : random_start + (fragment_size*5)]



        print "Uncleaned alignment length: %i" %(len(per_chr_alignment[0]))
            # create alignment and reduce scaffold length to shortest scaffold
        if options.clean:
            print "cleaning alignment"
            processed_alignment = clean_scaffold(per_chr_alignment) # perform clean of aligned scaffold
            if not processed_alignment:
                print "Something went wrong while cleaning up alignment of scaffold %s" %record
                continue
        else:
            processed_alignment = per_chr_alignment

        processed_alignment = processed_alignment[:, 0:fragment_size]

        if (len(processed_alignment[0])) != fragment_size:
            print "Cleaned alignment length not equal to %i, is %i" %(fragment_size, len(processed_alignment[0]))
            print "skipping..."
            continue

        with open(processed_alignment[0].id + "_" + str(fragment_size) + "_" + str(fragment_count) + ".fasta", "w") as fout:
            print "writing %s..." % (processed_alignment[0].id + str(fragment_size) + str(fragment_count))
            AlignIO.write(processed_alignment, fout, "fasta")
            print "done."

        fragment_count += 1
    return seqs


def main():
    global f_ab, f_extra, bt_positions
    seqs = {}

    records = []
    fname_list = [basename(fpath) for fpath in options.input_files]
    query_list = load_BED(options.bed_file) if options.bed_file else []

    assert type(options.fragment_size == "int")



    for fpath in options.input_files:  # read input fasta files
        print "loading fasta files..."
        fname = basename(fpath)
        seqs[fname] = SeqIO.index(fpath, "fasta")  # load sequences of each fasta file as SeqIO.index
        records_per_fasta = seqs.get(fname).keys()  # create list of fasta headers (i.e. scaffolds)
        records.extend(records_per_fasta)
        print "..."+fname

    print "\n"

    records = set(records)  # make unique list of scaffolds, i.e. remove duplicates

    #bt_positions = BedToolPositions()
    if options.random:
        for step in xrange(options.min, options.max, options.fragment_size):
            fragments = iterate_random(records, fname_list, query_list, seqs, step)
    elif not options.random:
        fragments = iterate_all(records, fname_list, query_list, seqs)




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
    parser.add_argument('-b', '--bed_file', required=False, help='Query file')
    parser.add_argument('-c', '--clean', required=False, action='store', help='Set if scaffolds should be cleaned')
    parser.add_argument('-r', '--random', required=False, type=int, help='Get X random fragments from the genome. ', default=0)
    parser.add_argument('--min', required=False, type=int, help='Minimum fragment size', default=0)
    parser.add_argument('--max', required=False, type=int, help='Maximum fragment size', default=10000)

    options = parser.parse_args()


    main()
