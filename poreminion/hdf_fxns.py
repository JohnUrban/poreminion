import sys
import re
import h5py
import seaborn
import numpy as np
from glob import glob
from Bio import SeqIO
##from cStringIO import StringIO
from collections import defaultdict
try:
    from cStringIO import *
except Exception:
    from StringIO import *


######################
''' Begin params ''' 
workflow = 'Basecall_2D_000'  # The workflow that was ran
metrichor_folder = 'data'# The directory containing the fast5 files from Metrichor
# #### Iterate over all the fast5 files and get the quality scores
files = glob('{}/*.fast5'.format(metrichor_folder))
''' End Params '''
######################



def analysesTypes(files):
    analyses = set()
    for fname in files:

        fh = h5py.File(fname, 'r')    
        
        key = 'Analyses'   
        if not key in fh:
            fh.close()
            continue
        else:
            for e in fh[key]:
                analyses.add(e)
        fh.close()
        
    return analyses



def get_events(f5,f5type="minknow",metrichor=None):
    if f5type == "metrichor":
        path = "Analyses/EventDetection_000/Reads/Read_0/Events"
    elif f5type == "minknow":
        path = "Analyses/EventDetection_000/Reads/Read_0/Events"
    for event in f5[path]:
        sys.stdout(event+"\n")


def get_num_events(f5,f5type="minknow"):
    if f5type == "metrichor":
        path = "Analyses/EventDetection_000/Reads/Read_0/Events"
    elif f5type == "minknow":
        path = "Analyses/EventDetection_000/Reads/Read_0/Events"
##    for event in f5[path]:
##        sys.stdout(event+"\n")
    sys.stdout(f5[path].len()+"\n")



def numEvents_vs_seqLen(f5):
    lengths = []
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    try: #input_event_len
        lengths.append(str(f5["Analyses/EventDetection_000/Reads/"+read+"/Events"].len()))
    except KeyError:
        lengths.append("-")
    try: # 2D_fastq_path
        lengths.append(str(len(str(SeqIO.read(StringIO(f5["Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"][()]), "fastq").seq))))
    except KeyError:
        lengths.append("-")
    try: ## 2D ratio
        lengths.append(str(float(lengths[-2])/float(lengths[-1])))
    except ValueError:
        lengths.append("-")
    try: # template_event_len
        lengths.append(str(f5["Analyses/Basecall_2D_000/BaseCalled_template/Events"].len()))
    except KeyError:
        lengths.append("-")
    try: # template_fastq_len 
        lengths.append(str(len(str(SeqIO.read(StringIO(f5["Analyses/Basecall_2D_000/BaseCalled_template/Fastq"][()]), "fastq").seq))))
    except KeyError:
        lengths.append("-")
    try: ## template ratio
        lengths.append(str(float(lengths[-2])/float(lengths[-1])))
    except ValueError:
        lengths.append("-")
    try: # complement_event_len
        lengths.append(str(f5["Analyses/Basecall_2D_000/BaseCalled_complement/Events"].len()))
    except KeyError:
        lengths.append("-")
    try: # complement_fastq_len
        lengths.append(str(len(str(SeqIO.read(StringIO(f5["Analyses/Basecall_2D_000/BaseCalled_complement/Fastq"][()]), "fastq").seq))))
    except KeyError:
        lengths.append("-")
    try: ## complement ratio
        lengths.append(str(float(lengths[-2])/float(lengths[-1])))
    except ValueError:
        lengths.append("-")
##    sys.stdout(("\t").join(lengths)+"\n")
    print ("\t").join(lengths)


def numEvents_vs_seqLen_iter(f5dir):
    files = glob(f5dir + "/" + "*fast5")
    for f in files:
        f5 = h5py.File(f)
        numEvents_vs_seqLen(f5)
        f5.close()
        
def printname(name):
    ''' f.visit(printname) will iterate and print all group/paths in hdf5 file where f is a hdf file object'''
    sys.stdout.write(name+'\n')

def groupattr(name):
    try:
        sys.stdout.write(str(name)+"\t"+str(name.attrs.items())+"\n")
    except AttributeError:
        sys.stdout.write(name+"\tNo Attrs\n")
        
def exploreHDFFileStructure(files):
    '''files is a list of filenames created by: glob('{}/*.fast5'.format(metrichor_folder))'''
    global groups
    groups = set()
    
    for fname in files:
        fh = h5py.File(fname, 'r')    
        fh.visit(lambda group: groups.add(group))
        fh.close()
    for s in sorted(list(groups)):
        sys.stdout.write(s+'\n')
        

def showHDFelement(files, group):
    for fname in files:
        fh = h5py.File(fname, 'r')
        if group in fh:
            try:
                for e in fh[group]:
                    print e
            except (TypeError, IOError):
                pass ##TOOD: handle each error differently
        else:
            sys.stderr.write("File does not contain this group.\n")
        fh.close()

def exploreHDFfileGroupAttrs(files):
    if len(files) > 1:
        print "For now just use this on a single file still in a list though [file]..."
    for fname in files:
        fh = h5py.File(fname, 'r')
        sys.stdout.write(str(fh.filename)+"\t"+str(fh.attrs.items())+"\n")
        fh.visit(groupattr)
        fh.close()

