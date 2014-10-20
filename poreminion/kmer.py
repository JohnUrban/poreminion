import sys
from poretools.Fast5File import *
import numpy as np
import pandas
import matplotlib.pyplot as plt
from collections import defaultdict
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

## attempt 1: brute force
def kmercount(string, kmerdict, k):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    for i in range(len(string)-k+1):
        kmerdict[string[i:i+k]] += 1
    return kmerdict

def writekmer(kmerdict, fileobj=sys.stdout):
    ''' kmerdict is a dict/default dict object'''
    for kmer in sorted(kmerdict.keys()): ## this could be a big memory suck
        fileobj.write(kmer + '\t' + str(kmerdict[kmer]) + '\n')


def run(parser, args):
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
                    print seqLen, fa.seq[:10] ## DELETE this line
                    kmerdict = kmercount(fa.seq, kmerdict, args.k)
            files_processed += 1
        if files_processed % 100 == 0:
                logger.info("%d files processed." % files_processed)
        fast5.close()
    if args.saveas:
        fhout = open(args.saveas, 'w')
        writekmer(kmerdict, fhout)
        fhout.close()
    else:
        writekmer(kmerdict)

        
    
