import sys
from poretools.Fast5File import *
from events_tools_pm import *
import numpy as np
from findTimeErrors import has_time_error, has_time_error_list


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
    '''readtpye in 2d, template, complement'''
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



def get_frag_stats(fast5):
        # name, numinputevents, fragsize-est, 2dreadlen, templen, complen, meanqual,
        f=fast5.hdf5file
        name = fast5.filename.split("/")[-1]
        numevents = get_num_events(f, "input")
        if fast5.has_2D():
            has2d=True
            seqlen2d = get_seq_len(f, "2d")
            fragsize = seqlen2d
            meanscore2d = get_mean_qscore(f, "2d")
        else:
            has2d=False
            seqlen2d = "-"
            meanscore2d = "-"
        seqlentemp = get_seq_len(f, "template")
        numtempevents = get_num_events(f, "template")
        meanscoretemp = get_mean_qscore(f, "template")
        if has_complement(f):
            hascomp=True
            seqlencomp = get_seq_len(f, "complement")
            numcompevents = get_num_events(f, "complement")
            meanscorecomp = get_mean_qscore(f, "complement")
            if not fast5.has_2D():
                if seqlencomp > seqlentemp:
                    fragsize = seqlencomp
                else:
                    fragsize = seqlentemp
        else:
            hascomp=False
            seqlencomp = "-"
            numcompevents = "-"
            meanscorecomp = "-"
            if not fast5.has_2D():
                fragsize = seqlentemp
        ## tc_ratio
        if hascomp:
            tc_ratio = str(float(numtempevents)/float(numcompevents))
        else:
            tc_ratio = "-"
        ## is fragsize estimate robust to all seq lengths?
        upper = fragsize+fragsize*0.2
        lower = fragsize-fragsize*0.2
        score = 0
        if seqlentemp <= upper and seqlentemp >= lower:
            score += 1
        if hascomp:
            if seqlencomp <= upper and seqlencomp >= lower:
                score += 1
        else: ## if no comp, then only temp, and fragsize robust by default
            score += 1 ## no need to check seqlen2d since if present it is fragsize and if not, no penalty
        if score == 2:
            robust = 1
        else:
            robust = 0
        return [name, fragsize, numevents, hascomp, has2d, numtempevents, numcompevents, seqlen2d, seqlentemp, seqlencomp, meanscore2d, meanscoretemp, meanscorecomp, tc_ratio, robust]


## get time error input (opt)
## input mean, sd, median, min, max durations...
## get numtempevents, numtempevents_0move, numtempevents_1move, numtempevents_2move, numtempevents_3move,numtempevents_4move,numtempevents_5move,
##      mean duration, sd duration
## same for complement
## overall events slope  (event_n-event1)/(time_n-time_1)
## template and complement slopes
##        
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
    

def run(parser, args):
    for fast5 in Fast5FileSet(args.files):
        fragstats = [str(e) for e in get_frag_stats(fast5)]
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
        fast5.close()


