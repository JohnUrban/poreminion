import sys
import h5py
import numpy as np

# raw files:            start, length, mean, var
# base called files:    mean, stddev, start, length
# order of parameters differs in raw vs. basecalled files
# also, std dev used in 1, and variance used in other
# need to interpret which filetype reading and print out in basecalled format


def is_basecalled(f5):
    try:
        f5['/Analyses/Basecall_2D_000']
        return True
    except KeyError:
        return False

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
                

