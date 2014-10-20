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
    r = robjects.r
    r.library("ggplot2")
    grdevices = importr('grDevices')
    kmerdict = kmercount_in_table(args.table1)
    data = defaultdict(list)
    numKmers = len(kmerdict)
    for k in sorted(kmerdict.keys()):
        data['kmers'].append(k)
        data['counts'].append(kmerdict[k])
    plt.bar(left=range(1,numKmers+1), height=data['counts'], width=1.0)
    plt.show()
    
    
def run(parser, args):
    if args.mpl:
        singleTablePlot_mpl(parser, args)
    elif args.gg:
        singleTablePlot_gg(parser, args)

        
    
