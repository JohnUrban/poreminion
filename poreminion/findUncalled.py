import sys
from poretools.Fast5File import *
from glob import glob
import h5py
from collections import defaultdict
from subprocess import call
import os

def update_dict(d, numevents):
    if numevents < d["min"]:
        d["min"] = numevents
    elif numevents > d["max"]:
        d["max"] = numevents
    d["sum"] += numevents
    d["num"] += 1
    return d

def print_dict(d, msg):
    try:
        mean = float(d["sum"])/d["num"]
    except ZeroDivisionError:
        mean = (float(d["sum"]))/(d["num"]+1e-323)
    return msg, "numfiles:"+str(d["num"]), "minNumEvents:"+str(d["min"]), "maxNumEvents:"+str(d["max"]), "meanNumEvents:"+str(mean)

def write_out(filename, lines):
    out = open(filename,'w')
    out.writelines(("\n").join(lines)+"\n")
    out.close()

def write_stats(filename, notemp, toofewevents, toomanyevents, ioerror):
    out = open(filename,'w')
    out.write(("\n").join(print_dict(notemp, "No template data found.")) +"\n\n")
    out.write(("\n").join(print_dict(toofewevents, "Number of events is less than the allowed mimimum.")) +"\n\n")
    out.write(("\n").join(print_dict(toomanyevents, "Number of events is more than the allowed maximum.")) +"\n\n")
    out.write(("\n").join(["IOerror", "numfiles:"+str(len(ioerror))]) +"\n\n")
    out.write("All too-big-to-base-call sizes:\n")
    out.write(("\n").join([str(e) for e in sorted(toomanyevents["sizes"])])+"\n")
    out.close()


def runfail(cmd):
    if not 0 == os.system(cmd):
        sys.exit("Failed : %s " % cmd)

        
def move_files(filelist, dirname):
    print "DIRNAME: ", dirname
    if not os.path.exists(dirname):
##        print "Making dir...."
        runfail("mkdir "+dirname)
    for f in filelist:
##        print f
        newf = os.path.join(dirname, f.split("/")[-1])
        runfail("mv {} {}".format(f, newf))
            
    

def run(parser, args):
    nocall = defaultdict(list)
    total = 0
    numNotCalled = 0
    notemp = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toofewevents = {"min":float("inf"),"max":0,"sum":0,"num":0}
    toomanyevents = {"min":float("inf"),"max":0,"sum":0,"num":0, 'sizes':[]}
    for fast5 in Fast5FileSet(args.files):
        total += 1
        fas = fast5.get_fastas_dict()
        try:
            f=h5py.File(fast5.filename)
        except IOError:
            nocall["IOError"].append(fast5.filename)
            continue
        if len(fas) == 0:
            numNotCalled += 1
            log=f["/Analyses/Basecall_2D_000/Log"].value
            readnum = [e for e in f["/Analyses/EventDetection_000/Reads/"]][0]
            numevents = f["/Analyses/EventDetection_000/Reads/"+readnum+"/Events"].shape[0]
##            print fast5.filename, str(log.split("\n"))
            if "No template data found." in log:
                nocall["No template data found."].append(fast5.filename)
                notemp = update_dict(notemp, numevents)
            elif "Number of events is less than the allowed mimimum." in log:
                nocall["Number of events is less than the allowed mimimum."].append(fast5.filename)
                toofewevents = update_dict(toofewevents, numevents)
            elif "Number of events is more than the allowed maximum." in log:
                nocall["Number of events is more than the allowed maximum."].append(fast5.filename)
                toomanevents = update_dict(toomanyevents, numevents)
                toomanyevents["sizes"].append(numevents)
        f.close()
        fast5.close()

    if args.move:
##        print args.files[0]
        base = args.files[0].split("/")
        if base[-1] == "":
            base = base[:-2]
        else:
            base = base[:-1]
        base = ("/").join(base)
##        print base
        if nocall["No template data found."]:
            move_files(nocall["No template data found."], os.path.join(base, "notemplate"))
        if nocall["Number of events is less than the allowed mimimum."]:
            move_files(nocall["Number of events is less than the allowed mimimum."], os.path.join(base, "toofewevents"))
        if nocall["Number of events is more than the allowed maximum."]:
            move_files(nocall["Number of events is more than the allowed maximum."], os.path.join(base, "toomanyevents"))
        if nocall["IOError"]:
            move_files(nocall["IOError"], os.path.join(base, "IOError"))
            

    write_out(args.outprefix+".notemplate.txt", nocall["No template data found."])
    write_out(args.outprefix+".toofew.txt", nocall["Number of events is less than the allowed mimimum."])
    write_out(args.outprefix+".toomany.txt", nocall["Number of events is more than the allowed maximum."])
    write_out(args.outprefix+".IOerror.txt", nocall["IOError"])
    write_stats(args.outprefix+".stats.txt", notemp, toofewevents, toomanyevents, nocall["IOerror"])
        
