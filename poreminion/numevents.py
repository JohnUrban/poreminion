import sys
from poretools.Fast5File import *
import h5py

            
    

def run(parser, args):
    notemp = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toofewevents = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toomanyevents = {"min":float("inf"),"max":0,"sum":0,"num":0, 'sizes':[]}
    for fast5 in Fast5FileSet(args.files):
        f=h5py.File(fast5.filename)
        readnum = [e for e in f["/Analyses/EventDetection_000/Reads/"]][0]
        numevents = f["/Analyses/EventDetection_000/Reads/"+readnum+"/Events"].shape[0]
        f.close()
        print fast5.filename + "\t" + str(numevents)
        fast5.close()

