import sys
from poretools.Fast5File import *
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
from collections import defaultdict
#
from kmer import kmercount_in_table

## attempt 1: brute force
def singleTablePlot_gg(parser, args):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    r = robjects.r
    r.library("ggplot2")
    grdevices = importr('grDevices')
    kmerdict = kmercount_in_table(args.table1)
    data = defaultdict(list)
    numKmers = len(kmerdict)
    for k in sorted(kmerdict.keys()):
        data['kmers'].append(k)
        data['counts'].append(kmerdict[k])
    df = robjects.DataFrame(data)
    gp = ggplot2.ggplot(df)
##    pp = gp + ggplot2.geom_bar(stat="identity")
    pp = gp + ggplot2.aes_string(x=range(1,numKmers+1),y=data['counts']) \
         + ggplot2.geom_bar(stat="identity") \
         + ggplot2.scale_x_continuous(name="kmer", breaks=0.5+(range(1,numKmers+1)), labels=kmers)
    pp.plot()
    print('Type enter to exit.')
    raw_input()
    
def singleTablePlot_mpl(parser, args):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    kmerdict = kmercount_in_table(args.table1)
    data = defaultdict(list)
    numKmers = len(kmerdict)
    for k in sorted(kmerdict.keys()):
        data['kmers'].append(k)
        data['counts'].append(kmerdict[k])
    plt.bar(left=range(1,numKmers+1), height=data['counts'], width=1.0)
    plt.show()
    
def ensureEqualKmerSets(kmerdict1, kmerdict2):
    ''' dicts are defaultdict(int)'''
    ## make sure they have same set of kmers, for any missing elements in set A or B, add it with count of 0
    for key in kmerdict1:
        kmerdict2[key]
    for key in kmerdict2:
        kmerdict1[key]
    return kmerdict1, kmerdict2


def kmerDictToPlotData(kmerdict):
    data = defaultdict(list)
    for k in sorted(kmerdict.keys()):
        data['kmers'].append(k)
        data['counts'].append(kmerdict[k])
    return data

def totalKmerCount(kmerdict):
    total = 0
    for v in kmerdict.values():
        total += v
    return total

def readInTwoKmerTables(parser, args):
    kmerdict1 = kmercount_in_table(args.table1)
    kmerdict2 = kmercount_in_table(args.table2)
    ## ensure equal kmer sets
    kmerdict1, kmerdict2 = ensureEqualKmerSets(kmerdict1, kmerdict2)
    return kmerdict1, kmerdict2

def twoTablePlot_mpl(parser, args):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    ## read in and equilibrate the 2 kmer count tables
    kmerdict1, kmerdict2 = readInTwoKmerTables(parser, args)
    
    ## make approp data structures
    data1 = kmerDictToPlotData(kmerdict1)
    data2 = kmerDictToPlotData(kmerdict2)

    plt.scatter(x=data1['counts'], y=data2['counts'], s=10)
    for i in range(len(data1['kmers'])):
        plt.annotate(data1['kmers'][i], (data1['counts'][i],data2['counts'][i]), size='xx-small')
    if args.saveas is not None and (args.saveas.endswith(".pdf") or args.saveas.endswith(".jpg")):
            plt.savefig(args.saveas)
    else:
        plt.show()

    
def run(parser, args):
    if args.table1 and not args.table2:
        if args.mpl:
            singleTablePlot_mpl(parser, args)
        elif args.gg:
            singleTablePlot_gg(parser, args)
    elif args.table1 and args.table2:
        twoTablePlot_mpl(parser, args)
        

        
    
