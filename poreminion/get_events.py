import h5py
from events_tools_pm import *
from fragstats import has_complement

def run(parser, args):
    f5 = h5py.File(args.fast5)
    if args.raw:
        basecalled = False
    elif args.basecalled:
        basecalled = True
    else:
        basecalled = is_basecalled(f5)

    if args.type == 'input':
        if args.header:
            print ("\t").join(['#mean', 'stddev', 'start', 'length'])
        get_input_events(f5, basecalled)
    elif args.type == 'template':
        if args.header and basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length', 'model_state', 'model_level', 'move', 'p_model_state', 'mp_state', 'p_mp_state', 'p_A', 'p_C', 'p_G', 'p_T'])
        elif args.header and not basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length'])
        get_template_events(f5, basecalled)
    elif args.type == 'complement':
        if has_complement(f5):
            if args.header and basecalled:
                print ("\t").join(['#mean', 'stddev', 'start', 'length', 'model_state', 'model_level', 'move', 'p_model_state', 'mp_state', 'p_mp_state', 'p_A', 'p_C', 'p_G', 'p_T'])
            elif args.header and not basecalled:
                print ("\t").join(['#mean', 'stddev', 'start', 'length'])
            get_complement_events(f5, basecalled)
        else:
            sys.stderr.write("Poreminion: No complement found.\n")
