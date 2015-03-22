import sys
from poretools.Fast5File import *
from fragstats import get_frag_stats
import h5py
from findUncalled import runfail, write_out, move_files, get_base_name
from subprocess import call
import os
#logging
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)


def has_time_error(f5connection, basecalled):
    ''' times is a list of numbers'''
    read = [e for e in f5connection["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    if basecalled:
        i=2
    else:
        i=0
    t1=f5connection[path][0][i]
    for event in f5connection[path][1:]:
        t2 = event[i]
        if t2 < t1:
            return True
        t1 = t2
    return False

def is_basecalled(f5connection):
    try:
        f5connection['/Analyses/Basecall_2D_000']
        return True
    except KeyError:
        return False

def run(parser, args):
    errorfiles = []
    fileset = Fast5FileSet(args.files)
    numfiles = float(fileset.num_files_in_set)
    files_processed = 0
    proportion_files = 0.1
    for fast5 in fileset:
        if has_time_error(fast5.hdf5file, is_basecalled(fast5.hdf5file)):
            errorfiles.append(fast5.filename)
            if not args.outprefix:
                print fast5.filename.split("/")[-1]
        fast5.close()
        files_processed += 1.0
        if files_processed/numfiles >= proportion_files and args.verbose:
                logger.info("%d percent files processed." % (100*files_processed/numfiles))
                proportion_files += 0.1
    if args.outprefix:
        write_out(args.outprefix+".time.errors.txt", errorfiles)
    if args.move:
        base = get_base_name(args.files[0])
        move_files(errorfiles, os.path.join(base, "time_errors"))
