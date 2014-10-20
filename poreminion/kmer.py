import sys
from poretools.Fast5File import *
import numpy as np
import pandas
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO

## attempt 1: brute force
def kmercount_in_string(string, kmerdict, k):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    for i in range(len(string)-k+1):
        kmerdict[string[i:i+k]] += 1
    return kmerdict

def writekmer(kmerdict, fileobj=sys.stdout):
    ''' kmerdict is a dict/default dict object'''
    total = float(sum(kmerdict.values()))
    for kmer in sorted(kmerdict.keys()): ## this could be a big memory suck
        fileobj.write(kmer + '\t' + str(kmerdict[kmer]) + '\t' + str(kmerdict[kmer]/total) + '\n')


def kmercount_in_fast5(parser, args):
    files_processed = 0
    kmerdict = defaultdict(int)
    for fast5 in Fast5FileSet(args.files):
        if args.start_time or args.end_time:
                read_start_time = fast5.get_start_time()
                read_end_time = fast5.get_end_time()
                if args.start_time and args.start_time > read_start_time:
                        fast5.close()
                        continue
                if args.end_time and args.end_time < read_end_time:
                        fast5.close()
                        continue
        if args.high_quality:
            if fast5.get_complement_events_count() <= \
               fast5.get_template_events_count():
                    fast5.close()
                    continue
        if args.single_read:
            fas = [fast5.get_fasta()]
        else:
            fas = fast5.get_fastas(args.type)
        for fa in fas:
            seqLen = len(fa.seq)
            if fa is not None and not (seqLen < args.min_length or seqLen > args.max_length):
                    kmerdict = kmercount_in_string(fa.seq, kmerdict, args.k)
            files_processed += 1
        if files_processed % 100 == 0:
                logger.info("%d files processed." % files_processed)
        fast5.close()
    return kmerdict


def kmercount_in_fastx(parser, args, fastx="fasta"):
    kmerdict = defaultdict(int)
    if fastx == "fasta":
        fh = args.fasta
    elif fastx == "fastq":
        fh = args.fastq
    for fa in SeqIO.parse(fh, fastx):
        seqLen = len(fa.seq)
        if fa is not None and not (seqLen < args.min_length or seqLen > args.max_length):
                kmerdict = kmercount_in_string(str(fa.seq), kmerdict, args.k)
                if args.rev_comp:
                    kmerdict = kmercount_in_string(str(fa.reverse_complement().seq), kmerdict, args.k)
        
    return kmerdict


def kmercount_in_table(kmerCountFile):
    ''' Input is a kmer count table.'''
    kmerdict = defaultdict(int)
    fh = kmerCountFile
    for line in open(fh, 'r'):
        kmer, count, prop = line.rstrip().split('\t')
        kmerdict[kmer] = int(count)
    return kmerdict

def run(parser, args):
    if args.fasta:
        kmerdict = kmercount_in_fastx(parser, args, fastx="fasta")
    elif args.fastq:
        kmerdict = kmercount_in_fastx(parser, args, fastx="fastq")
    else:
        kmerdict = kmercount_in_fast5(parser, args)
    if args.saveas:
        fhout = open(args.saveas, 'w')
        writekmer(kmerdict, fhout)
        fhout.close()
    else:
        writekmer(kmerdict)

        
    
