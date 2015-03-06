import sys
from poretools.Fast5File import *
import h5py

            
def get_num_events(f5connection):
        readnum = [e for e in f5connection["/Analyses/EventDetection_000/Reads/"]][0]
        numevents = f5connection["/Analyses/EventDetection_000/Reads/"+readnum+"/Events"].shape[0]
        return numevents

def run(parser, args):
    notemp = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toofewevents = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toomanyevents = {"min":float("inf"),"max":0,"sum":0,"num":0, 'sizes':[]}
    for fast5 in Fast5FileSet(args.files):
        numevents = get_num_events(fast5.hdf5file)
        print fast5.filename + "\t" + str(numevents)
        fast5.close()

