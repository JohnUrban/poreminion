import sys
from poretools.Fast5File import *

def run(parser, args):
    total = 0
    twodir = 0
    for fast5 in Fast5FileSet(args.files):
        if args.only2d:
            if fast5.has_2D():
                print fast5.filename
        elif args.no2d:
            if not fast5.has_2D():
                print fast5.filename
        else:
            print fast5.filename + "\t" + str(fast5.has_2D())
        fast5.close()


