import sys
from poretools.Fast5File import *
from events_tools_pm import *
import numpy as np
from findTimeErrors import has_time_error, has_time_error_list
import pandas
from nX import NX
from joblib import Parallel, delayed ## parallelize 0.3.2
#logging
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)

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

def has_complement(f5connection):
    try:
        f5connection["/Analyses/Basecall_2D_000/BaseCalled_complement"]
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

def get_basename(f5):
    return f5["/Analyses/Basecall_2D_000/Configuration/general"].attrs["basename"]

def get_channel_number(f5):
    return f5['/UniqueGlobalKey/channel_id'].attrs['channel_number']

def get_heatsink_temp(f5):
    return f5['/UniqueGlobalKey/tracking_id'].attrs['heatsink_temp'] 


def get_frag_stats(fast5):
        # name, numinputevents, fragsize-est, 2dreadlen, templen, complen, meanqual,
        f=fast5.hdf5file
        name = fast5.filename.split("/")[-1]
        numevents = get_num_events(f, "input")
        channelnum = get_channel_number(f)
        heatsink = get_heatsink_temp(f)
        if fast5.has_2D():
            has2d=1
            seqlen2d = get_seq_len(f, "2d")
            fragsize = seqlen2d
            meanscore2d = get_mean_qscore(f, "2d")
        else:
            has2d=0
            seqlen2d = "-"
            meanscore2d = "-"
        seqlentemp = get_seq_len(f, "template")
        numtempevents = get_num_events(f, "template")
        meanscoretemp = get_mean_qscore(f, "template")
        numcalledeventstemp = get_num_called_events(f, "template")
        numskipstemp = get_num_skips(f, "template")
        numstaystemp = get_num_stays(f, "template")
        strandscoretemp = get_strand_score(f, "template")
        numstutterstemp = get_num_stutters(f, "template")
        if has_complement(f):
            hascomp=1
            seqlencomp = get_seq_len(f, "complement")
            numcompevents = get_num_events(f, "complement")
            meanscorecomp = get_mean_qscore(f, "complement")
            numcalledeventscomp = get_num_called_events(f, "complement")
            numskipscomp = get_num_skips(f, "complement")
            numstayscomp = get_num_stays(f, "complement")
            strandscorecomp = get_strand_score(f, "complement")
            numstutterscomp = get_num_stutters(f, "complement")
            if not fast5.has_2D():
                if seqlencomp > seqlentemp:
                    fragsize = seqlencomp
                else:
                    fragsize = seqlentemp
        else:
            hascomp=0
            seqlencomp = "-"
            numcompevents = "-"
            meanscorecomp = "-"
            numcalledeventscomp = "-"
            numskipscomp = "-"
            numstayscomp = "-"
            strandscorecomp = "-"
            numstutterscomp = "-"
            if not fast5.has_2D():
                fragsize = seqlentemp
        ## tc_ratio --- QUESTION: which is more appropriate here: numcalledevents or numevents? Either way, both can be calculated afterward.
        if hascomp:
            log2_tc_ratio = str(np.log2(float(numtempevents)/float(numcompevents)))
        else:
            log2_tc_ratio = "-" 
        return [name, fragsize, numevents, hascomp, has2d, numtempevents, numcompevents, seqlen2d, seqlentemp, seqlencomp, meanscore2d, meanscoretemp, meanscorecomp, log2_tc_ratio, channelnum, heatsink, numcalledeventstemp, numcalledeventscomp, numskipstemp, numskipscomp, numstaystemp, numstayscomp, strandscoretemp, strandscorecomp, numstutterstemp, numstutterscomp]


       
def get_more_fragstats(f5, hascomp, check_time = False):
    ## f5 is the connection to fast5 file
    basecalled = is_basecalled(f5)
    input_events = store_input_events(f5, basecalled)
    if check_time:
        time_error = [int(has_time_error_list(input_events['start']))]
    else:
        time_error = []
    start_time = input_events['start'][0]
    end_time = input_events['start'][-1]
    nevents = len(input_events['start'])
    slope = (end_time-start_time)/(nevents-1)
    mean_dur_input = np.mean(input_events['length'])
    sd_dur_input = np.std(input_events['length'])
    med_dur_input = np.median(input_events['length'])
    min_dur_input = min(input_events['length'])
    max_dur_input = max(input_events['length'])
    template_events = store_template_events(f5, basecalled, getbasecallinfo=True)
    tmoves = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
    for move in template_events['move']:
        tmoves[move] += 1
    if hascomp:
        comp_events = store_complement_events(f5, basecalled, getbasecallinfo=True)
        cmoves = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
        for move in comp_events['move']:
            cmoves[move] += 1
    else:
        cmoves = {0:"-", 1:"-", 2:"-", 3:"-", 4:"-", 5:"-"}
    return [start_time, end_time, slope, mean_dur_input, sd_dur_input, med_dur_input, min_dur_input, max_dur_input, tmoves[0], tmoves[1], tmoves[2], tmoves[3], tmoves[4], tmoves[5], cmoves[0], cmoves[1], cmoves[2], cmoves[3], cmoves[4], cmoves[5]] + time_error
    


def get_empty_fragstats_dict(extensive=False, timecheck=False):
    fragstats = {'name':[], 'fragsize':[], 'numevents':[], 'hascomp':[], 'has2d':[], 'numtempevents':[], 'numcompevents':[], 'seqlen2d':[], 'seqlentemp':[], 'seqlencomp':[], 'meanscore2d':[], 'meanscoretemp':[], 'meanscorecomp':[], 'log2_tc_ratio':[], 'channelnum':[], 'heatsink':[], 'numcalledeventstemp':[], 'numcalledeventscomp':[], 'numskipstemp':[], 'numskipscomp':[], 'numstaystemp':[], 'numstayscomp':[], 'strandscoretemp':[], 'strandscorecomp':[], 'numstutterstemp':[],'numstutterscomp':[]}
    if extensive:
        fragstats.update({'starttime':[], 'endtime':[], 'slope':[], 'mean_dur_input':[], 'sd_dur_input':[], 'med_dur_input':[], 'min_dur_input':[], 'max_dur_input':[], 'tmove_0':[], 'tmove_1':[], 'tmove_2':[], 'tmove_3':[],'tmove_4':[], 'tmove_5':[], 'cmove_0':[], 'cmove_1':[], 'cmove_2':[], 'cmove_3':[], 'cmove_4':[], 'cmove_5':[]})
    if timecheck:
        fragstats.update({'timeerror':[]})
    return fragstats

def update_fragstats_dict(fragstats, newfrag, extensive=False, timecheck=False):
    ## fragstats is a fragstats dict -- e.g. output of get_empty_fragstats_dict()
    ## newfrag is a line.split() from a fragstats tab-delim table file
    ## extensive is T/F whether or not the fragstats file is output of --extensive
    ## timecheck is T/F whether the fragstats file had the --timecheck used
    fragstats['name'].append(newfrag[0])
    fragstats['fragsize'].append(float(newfrag[1]))
    fragstats['numevents'].append(float(newfrag[2]))
    fragstats['hascomp'].append(int(newfrag[3]))
    fragstats['has2d'].append(int(newfrag[4]))
    fragstats['numtempevents'].append(float(newfrag[5]))
    fragstats['seqlentemp'].append(float(newfrag[8]))
    fragstats['meanscoretemp'].append(float(newfrag[11]))
    fragstats['numcalledeventstemp'].append(float(newfrag[16]))
    fragstats['numskipstemp'].append(float(newfrag[18]))
    fragstats['numstaystemp'].append(float(newfrag[20]))
    fragstats['strandscoretemp'].append(float(newfrag[22]))
    if newfrag[24] != "-": ##numstutters only in newer fast5s
        fragstats['numstutterstemp'].append(float(newfrag[24]))
    else:
        fragstats['numstutterstemp'].append(None)
    if fragstats['hascomp'][-1]:
        fragstats['numcompevents'].append(float(newfrag[6]))
        fragstats['seqlencomp'].append(float(newfrag[9]))
        fragstats['meanscorecomp'].append(float(newfrag[12]))
        fragstats['log2_tc_ratio'].append(float(newfrag[13]))
        fragstats['numcalledeventscomp'].append(float(newfrag[17]))
        fragstats['numskipscomp'].append(float(newfrag[19]))
        fragstats['numstayscomp'].append(float(newfrag[21]))
        fragstats['strandscorecomp'].append(float(newfrag[23]))
        if newfrag[25] != "-":
            fragstats['numstutterscomp'].append(float(newfrag[25]))
        else:
            fragstats['numstutterscomp'].append(None)
    else:
        fragstats['numcompevents'].append(None)
        fragstats['seqlencomp'].append(None)
        fragstats['meanscorecomp'].append(None)
        fragstats['log2_tc_ratio'].append(None)
        fragstats['numcalledeventscomp'].append(None)
        fragstats['numskipscomp'].append(None)
        fragstats['numstayscomp'].append(None)
        fragstats['strandscorecomp'].append(None)
        fragstats['numstutterscomp'].append(None)
    if fragstats['has2d'][-1]:
        fragstats['seqlen2d'].append(float(newfrag[7]))
        fragstats['meanscore2d'].append(float(newfrag[10]))
    else:
        fragstats['seqlen2d'].append(None)
        fragstats['meanscore2d'].append(None)

    fragstats['channelnum'].append(float(newfrag[14]))
    fragstats['heatsink'].append(float(newfrag[15]))
    
    if extensive:
        fragstats['starttime'].append(float(newfrag[26]))
        fragstats['endtime'].append(float(newfrag[27]))
        fragstats['slope'].append(float(newfrag[28]))
        fragstats['mean_dur_input'].append(float(newfrag[29]))
        fragstats['sd_dur_input'].append(float(newfrag[30]))
        fragstats['med_dur_input'].append(float(newfrag[31]))
        fragstats['min_dur_input'].append(float(newfrag[32]))
        fragstats['max_dur_input'].append(float(newfrag[33]))
        fragstats['tmove_0'].append(float(newfrag[34]))
        fragstats['tmove_1'].append(float(newfrag[35]))
        fragstats['tmove_2'].append(float(newfrag[36]))
        fragstats['tmove_3'].append(float(newfrag[37]))
        fragstats['tmove_4'].append(float(newfrag[38]))
        fragstats['tmove_5'].append(float(newfrag[39]))
        if fragstats['hascomp'][-1]:
            fragstats['cmove_0'].append(float(newfrag[40]))
            fragstats['cmove_1'].append(float(newfrag[41]))
            fragstats['cmove_2'].append(float(newfrag[42]))
            fragstats['cmove_3'].append(float(newfrag[43]))
            fragstats['cmove_4'].append(float(newfrag[44]))
            fragstats['cmove_5'].append(float(newfrag[45]))
        else:
            fragstats['cmove_0'].append(None)
            fragstats['cmove_1'].append(None)
            fragstats['cmove_2'].append(None)
            fragstats['cmove_3'].append(None)
            fragstats['cmove_4'].append(None)
            fragstats['cmove_5'].append(None)
    if timecheck and not extensive:
        fragstats['timeerror'].append(float(newfrag[26]))
    elif timecheck and extensive:
        fragstats['timeerror'].append(float(newfrag[46]))
    return fragstats
    


def make_fragstats_dict(fragstatsfile, extensive=False, timecheck=False):
    ## fragstatsfile is a fragstats tab-delim table file
    fragstats = get_empty_fragstats_dict(extensive, timecheck)
    for line in open(fragstatsfile,'r'):
        newfrag = line.strip().split("\t")
        fragstats = update_fragstats_dict(fragstats, newfrag, extensive, timecheck)
    return fragstats

def make_fragstats_dataframe(fragstatsfile, extensive=False, timecheck=False):
    ## TODO(?): could actually use pandas.read_table()
    return pandas.DataFrame(make_fragstats_dict(fragstatsfile, extensive, timecheck))



def plot_fragstats():
    pass


def frag_fxn(fast5, args):
    fragstats = get_frag_stats(fast5)
    if args.extensive:
        hascomp = has_complement(fast5.hdf5file)
        if args.checktime:
            fragstats += get_more_fragstats(fast5.hdf5file, hascomp = hascomp, check_time = True)
        else:
            fragstats += get_more_fragstats(fast5.hdf5file, hascomp = hascomp, check_time = False)
    if not args.extensive and args.checktime:
        input_events = store_input_events(fast5.hdf5file, basecalled)
        fragstats += [int(has_time_error_list(input_events['start']))]
    print ("\t").join([str(e) for e in fragstats])


def run(parser, args):
    pass
##    Parallel(n_jobs=args.parallel)(
##    delayed(frag_fxn)(fast5, args) for fast5 in Fast5FileSet(args.files))

    for fast5 in Fast5FileSet(args.files):
        frag_fxn(fast5, args)
        fast5.close()



