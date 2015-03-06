import sys
from poretools.Fast5File import *



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
        f=h5py.File(fast5.filename)
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
        f.close()
        return [name, fragsize, numevents, hascomp, has2d, numtempevents, numcompevents, seqlen2d, seqlentemp, seqlencomp, meanscore2d, meanscoretemp, meanscorecomp]
        

def run(parser, args):
    for fast5 in Fast5FileSet(args.files):
        print ("\t").join([str(e) for e in get_frag_stats(fast5)])
        fast5.close()


