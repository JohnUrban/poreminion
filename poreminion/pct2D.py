import sys
from poretools.Fast5File import *

def run(parser, args):
    total = 0
    twodir = 0
    for fast5 in Fast5FileSet(args.files):
        total += 1
        twodir += fast5.has_2D()
        fast5.close()
        
    print ("\t").join([str(e) for e in [total, twodir, 100.0*twodir/total]])
