from info import *
from poretools import *


def run(parser, args):
    gettemp, getcomp, get2d = True, True, True
    if args.nottemp:
        gettemp=False
    if args.notcomp:
        getcomp=False
    if args.not2d:
        get2d=False
    for fast5 in Fast5FileSet(args.files):
        f5 = fast5.hdf5file
        print_seq_name_and_length(f5, gettemp=gettemp, getcomp=getcomp, get2d=get2d)
        fast5.close()
