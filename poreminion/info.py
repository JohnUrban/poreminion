#info
import h5py, os
import poreminion
#logging
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)

def is_basecalled(f5):
    try:
        f5['/Analyses/Basecall_2D_000']
        return True
    except KeyError:
        return False

def get_mean_qscore(f5connection, readtype):
    '''readtpye in 2d, template, complement'''
    if readtype == "2d":
        path="/Analyses/Basecall_2D_000/Summary/basecall_2d"
    elif readtype == "template":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_template"
    elif readtype == "complement":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_complement"
    return f5connection[path].attrs["mean_qscore"]

def get_seq_len(f5connection, readtype):
    '''readtpye in 2d, template, complement'''
    if readtype == "2d":
        path="/Analyses/Basecall_2D_000/Summary/basecall_2d"
    elif readtype == "template":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_template"
    elif readtype == "complement":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_complement"
    return f5connection[path].attrs["sequence_length"]


def print_seq_name_and_length(f5, gettemp=True, getcomp=True, get2d=True):
    ## seq name in format of g4 seq names as well as staypos seqnames
    ## assumes basecalled
    ## f5 is hdf5 file connection
    name = get_basename(f5)
    if name.endswith("_strand"):
        name = name[:-6]
    if gettemp:
        print name + "template\t" + str(get_seq_len(f5, readtype="template"))
    if has_complement(f5):
        if getcomp:
            print name + "complement\t" + str(get_seq_len(f5, readtype="complement"))
        if has_2d(f5) and get2d:
            print name + "2D\t" + str(get_seq_len(f5, readtype="2d"))

def has_complement(f5connection):
    try:
        f5connection["/Analyses/Basecall_2D_000/BaseCalled_complement"]
        return True
    except KeyError:
        return False

def has_2d(f5connection):
    try:
        f5connection["/Analyses/Basecall_2D_000/BaseCalled_2D"]
        return True
    except KeyError:
        return False

def get_num_events(f5connection, eventstype):
    '''Assumes exists'''
    '''eventstpye in 2d, template, complement'''
    if eventstype == "input":
        path="/Analyses/Basecall_2D_000/Summary/split_hairpin"
    elif eventstype == "template":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_template"
    elif eventstype == "complement":
        path="/Analyses/Basecall_2D_000/Summary/basecall_1d_complement"
    try:
        return f5connection[path].attrs["num_events"]
    except KeyError:
        return None

def get_num_called_events(f5, eventstype):
    '''Assumes exists'''
    '''eventstpye in template, complement (no input)'''
    if eventstype == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['called_events']
    elif eventstype == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['called_events']
    

def get_num_skips(f5, eventstype):
    '''Assumes exists'''
    '''eventstpye in template, complement (no input)'''
    if eventstype == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['num_skips']
    elif eventstype == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['num_skips']

def get_num_stays(f5, eventstype):
    '''Assumes exists'''
    '''eventstpye in template, complement (no input)'''
    if eventstype == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['num_stays']
    elif eventstype == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['num_stays']

def get_strand_score(f5, eventstype):
    '''Assumes exists'''
    '''eventstpye in template, complement (no input)'''
    if eventstype == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['strand_score']
    elif eventstype == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['strand_score']

def get_num_stutters(f5, eventstype):
    '''Assumes exists'''
    '''eventstpye in template, complement (no input)'''
    try:
        if eventstype == "template":
            return f5['/Analyses/Basecall_2D_000/Summary/post_process_template'].attrs['num_stutters_found']
        elif eventstype == "complement":
            return f5['/Analyses/Basecall_2D_000/Summary/post_process_complement'].attrs['num_stutters_found']
    except: ##numstutters only in newer fast5s
        logger.info(get_basename(f5) + " does not have the 'post_process_" + eventstype + " groups")
        return "-"

def get_channel_number(f5):
    return f5['/UniqueGlobalKey/channel_id'].attrs['channel_number']

def get_heatsink_temp(f5):
    return f5['/UniqueGlobalKey/tracking_id'].attrs['heatsink_temp'] 


def get_minknow_version(f5):
    return f5['/UniqueGlobalKey/tracking_id'].attrs["version_name"]

def get_basecaller_events_limits(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['min_events'], f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['max_events']

def get_min_max_ratio(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/hairpin_align'].attrs['min_ratio'], f5['/Analyses/Basecall_2D_000/Configuration/hairpin_align'].attrs['max_ratio']

def get_model_type(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['model_type']


def get_basename(f5):
    return f5["/Analyses/Basecall_2D_000/Configuration/general"].attrs["basename"]

def get_sequence(f5, seqtype):
    ''' f5: f5connection
        seqtype: template, complement, 2d'''
    pass
    

def run_get_attributes(f5file):
    script = os.path.join(poreminion.scripts.__path__[0], 'getAttributes.sh')
    cmd = script + " " + f5file
    os.system(cmd)

## TODO: add more get functions
## TODO: add option allowing user to request specific attributes

def run(parser, args):
    if args.basic:
        f5 = h5py.File(args.fast5)
        basecalled = is_basecalled(f5)
        print get_basename(f5)
        print get_minknow_version(f5)
        if basecalled:
            print get_basecaller_events_limits(f5)
            print get_min_max_ratio(f5)
            print get_model_type(f5)
    if args.all:
        run_get_attributes(args.fast5)
