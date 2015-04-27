from events_tools_pm import *
from info import *
from fragstats import *

def staypos_fxn(f5, args):
    ## assumes basecalled
    ## f5 is hdf5 file connection
    name = get_basename(f5)
    if name.endswith("_strand"):
        name = name[:-6]
    events = store_template_events(f5, basecalled=True, getbasecallinfo=True)
    stay_locations_BED(events, name+"template")
    if has_complement(f5):
        events = store_complement_events(f5, basecalled=True, getbasecallinfo=True)
        stay_locations_BED(events, name+"complement")


def run(parser, args):
    for fast5 in Fast5FileSet(args.files):
        f5 = fast5.hdf5file
        staypos_fxn(f5, args)
        fast5.close()    
