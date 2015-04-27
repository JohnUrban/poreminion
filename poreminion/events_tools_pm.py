import sys
import h5py
import numpy as np
from info import *

# raw files:            start, length, mean, var
# base called files:    mean, stddev, start, length
# order of parameters differs in raw vs. basecalled files
# also, std dev used in 1, and variance used in other
# need to interpret which filetype reading and print out in basecalled format




def get_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    if basecalled:
        for event in f5[path]:
            print ("\t").join([str(e) for e in event])
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            print ("\t").join([str(e) for e in event])



def yield_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    if basecalled:
        for event in f5[path]:
            yield ("\t").join([str(e) for e in event])
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            yield ("\t").join([str(e) for e in event])


def update_events(events, newevent, getbasecallinfo=False):
    #mean, stddev, start, length
    events['mean'].append(float(newevent[0]))
    events['stdev'].append(float(newevent[1]))
    events['start'].append(float(newevent[2]))
    events['length'].append(float(newevent[3]))
    if getbasecallinfo:
        #mean, stddev, start, length, model_state, model_level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T
        events['model_state'].append(str(newevent[4]))
        events['model_level'].append(float(newevent[5]))
        events['move'].append(int(newevent[6]))
        events['p_model_state'].append(float(newevent[7]))
        events['mp_state'].append(str(newevent[8]))
        events['p_mp_state'].append(float(newevent[9]))
        events['p_A'].append(float(newevent[10]))
        events['p_C'].append(float(newevent[11]))
        events['p_G'].append(float(newevent[12]))
        events['p_T'].append(float(newevent[13]))
    return events


def store_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    events = {'mean':[], 'stdev':[], 'start':[], 'length':[]}
    if basecalled:
        for event in f5[path]:
            events = update_events(events, event)
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            events = update_events(events, event)
    for key in events.keys():
        events[key] = np.array(events[key])
    return events
            

def get_stranded_events(f5, strand="template"):
    ## mean, start, stddev, length, model state, model level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T
    ## with 1 observation of support, the stranded times are time/5000.0 --- for whatever reason (goes for start time and duration)
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "/Analyses/Basecall_2D_000/BaseCalled_" + strand + "/Events"
    for event in f5[path]:
        event = list(event)
        event = [event[0]] + [event[2]] + [event[1]] + event[3:]
        print ("\t").join([str(e) for e in event])

def store_stranded_events(f5, strand="template", getbasecallinfo=True):
    ## mean, start, stddev, length, model state, model level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T
    ## with 1 observation of support, the stranded times are time/5000.0 --- for whatever reason (goes for start time and duration)
    #when getbasecallinfo=True: mean, stddev, start, length, model_state, model_level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T   
    # else just: mean, stddev, start, length
    events = {'mean':[], 'stdev':[], 'start':[], 'length':[], 'model_state':[], 'model_level':[], 'move':[], 'p_model_state':[], 'mp_state':[], 'p_mp_state':[], 'p_A':[], 'p_C':[], 'p_G':[], 'p_T':[]}
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "/Analyses/Basecall_2D_000/BaseCalled_" + strand + "/Events"
    for event in f5[path]:
        event = list(event)
        event = [event[0]] + [event[2]] + [event[1]] + event[3:]
        events = update_events(events, event, getbasecallinfo=getbasecallinfo)
    for key in events.keys():
        events[key] = np.array(events[key])
    return events

def get_template_events(f5, basecalled):
    if basecalled:
        get_stranded_events(f5, strand="template")
    else:
        pass            

def store_template_events(f5, basecalled, getbasecallinfo=True):
    if basecalled:
        return store_stranded_events(f5, strand="template", getbasecallinfo=getbasecallinfo)
    else:
        getbasecallinfo=False #override
        pass

def get_complement_events(f5, basecalled):
    if basecalled:
        get_stranded_events(f5, strand="complement")
    else:
        pass

def store_complement_events(f5, basecalled, getbasecallinfo=True):
    if basecalled:
        return store_stranded_events(f5, strand="complement", getbasecallinfo=getbasecallinfo)
    else:
        getbasecallinfo=False
        pass


def reconstruct_sequence_from_stranded_events(events):
    ## assumes events are stranded/base-called
    length = len(events['mean'])
    seq = events['model_state'][0]
    for i in range(1,length):
        if events['move'][i] > 0:
            print events['move'][i] 
            seq += events['model_state'][i][5-int(events['move'][i]):]
    return seq

def model_state_start_index_of_stranded_events(events):
    ## assumes events are stranded/base-called
    ## returns updated events dict with seqindex key:value pairs
    length = len(events['mean'])
    events['seqindex'] = [0]
    for i in range(1,length):
        events['seqindex'].append(events['seqindex'][-1]+int(events['move'][i]))
    return events


def stay_locations_BED(events, name):
    ## assumes events are stranded/base-called
    ## prints BED entries of locations of 5mers that have move=0
    ## it will print same location as many times as there are stays at it
    ## by default position 0 is move 0 but is not a "stay", so start at 1
    length = len(events['mean'])
    index = 0
    for i in range(1,length):
        if events['move'][i] == 0:
            print ("\t").join([name, str(index), str(index+5)])
        else:
            index += int(events['move'][i])
