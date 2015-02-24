import sys
from poretools.Fast5File import *

def run(parser, args):
    total = 0
    twodir = 0
    for fast5 in Fast5FileSet(args.files):
        total += 1
        fas = fast5.get_fastas_dict()
        if len(fas) > 0:
            basecalled_files += 1
        if fast5.has_2D():
            
##        print len(fast5.get_fastas("template"))
        fast5.close()
        
    print ("\t").join([str(e) for e in [total, twodir, 100.0*twodir/total]])
