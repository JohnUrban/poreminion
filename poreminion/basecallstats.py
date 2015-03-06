import sys
from poretools.Fast5File import *
from fragstats import get_frag_stats

##Fragstats + distribution of moves + whether all events occur correctly in time
## make another function that plots dist of durations
##      as well as duration length vs. event num

def get_base_call_stats(fast5):
    basecall_stats = get_frag_stats(fast5)
    pass



def run(parser, args):
    for fast5 in Fast5FileSet(args.files):
        print ("\t").join([str(e) for e in get_frag_stats(fast5)])
        fast5.close()
