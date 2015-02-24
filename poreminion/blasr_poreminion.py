from subprocess import call
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

## plots:
# 1 - pctID vs. readLen
# 2 - score vs. readLen
# 3 - qCov vs pctID
# 4 - 

def blasr_call(reads, reference, out, return_stats=True):
    blasr = call(["blasr", reads, reference, "-out", out, "-m", "5"]) ## This output is annoying do not want last 3 col
    o = open(out, 'r')
    f = open(out+".tmp", 'w')
    colnames = "qName qLength qStart qEnd qStrand qAlnLen qCov tName tLength tStart tEnd tStrand tAlnLen tCov score numMatch numMismatch numIns numDel mapQV percentIdentity"
    f.write("#"+colnames+"\n")
    colnames = colnames.split()
    blasrstats = {colname:[] for colname in colnames}
    for line in o:
        line = line.strip().split()
        qAlnLen = str(int(line[3])-int(line[2]))
        qCov = str(100.0*(int(line[3])-int(line[2]))/int(line[1]))
        tAlnLen = str(int(line[8])-int(line[7]))
        tCov = str(100.0*(int(line[8])-int(line[7]))/int(line[6]))
        accuracy = str(100*float(line[11])/sum([int(line[11]), int(line[12]), int(line[13]), int(line[14])]))
        stats = line[:5]+[qAlnLen, qCov]+line[5:10]+[tAlnLen, tCov]+line[10:16]+[accuracy]
        f.write(("\t").join(stats)+"\n")
        if return_stats:
            for i in range(len(stats)):
                blasrstats[colnames[i]].append(stats[i])
    o.close()
    f.close()
    call(["mv", out+".tmp", out])
    if return_stats:
        return blasrstats


def read_blasr_output(blasrfile):
    f = open(blasrfile, 'r')
    colnames = "qName qLength qStart qEnd qStrand qAlnLen qCov tName tLength tStart tEnd tStrand tAlnLen tCov score numMatch numMismatch numIns numDel mapQV percentIdentity".split()
    blasrstats = {colname:[] for colname in colnames}
    for line in f:
        if line[0] != "#":
            stats = line.strip().split()
            for i in range(len(stats)):
                    blasrstats[colnames[i]].append(stats[i])
    return blasrstats

def reference_coverage(blasrstats):
    cov = defaultdict(list)
    numreads = len(blasrstats["tStart"])
    for i in range(numreads):
        cov[blasrstats["tName"][i]] += range(int(blasrstats["tStart"][i]), int(blasrstats["tEnd"][i])+1)
    return cov


def run(parser, args):
    # defualts
    need_to_align = True
    args.return_stats = False
    if args.plot:
        args.return_stats = True
    if args.blasr_file:
        need_to_align = False
    reference = "/Users/johnurban/data/otherGenomes/Lambda/lambdaGenomeSequence.fa"
    reads = "/Users/johnurban/searchPaths/github/poreminion/winner.2D.fa"
    reads = "Lambdarun001.fa"
    out = "del.blasr"
    if need_to_align:
        blasrstats = blasr_call(reads, reference, out, return_stats=args.return_stats)
    else:
        blasrstats = read_blasr_output(out)

    x = reference_coverage(blasrstats)
    heights,bins = np.histogram(x["lambdaGenome"], bins=100)
    norm = float(len(blasrstats["tStart"]))
    heights = [height/norm for height in heights] 
    centers = (bins[:-1]+bins[1:])/2.0
    width = bins[1]-bins[0]
    plt.bar(centers, heights, align = 'center', width=width)
##    plt.hist(x["lambdaGenome"], bins=100)


##    print b
    plt.show()
                   
##    print blasrstats
##    print b
##    plt.scatter(blasrstats['qLength'], blasrstats['percentIdentity'])
##    plt.show()
##                
    





##BLASR -- https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
## default output -- qName tName qStrand tStrand score percentSimilarity tStart tEnd tLength qStart qEnd qLength nCells
# space-delim -- 11 fields
# blasr winner.fa genome.fa  
# blasr winner.fa genome.fa -m 1

##SAM output
# blasr winner.fa genome.fa  -sam 
#  (1) "XS": read alignment start position without counting previous soft clips 
#  (2) "XE": read alignment end position without counting previous soft clips 
#  (3) "XL": aligned read length
#  (4) "XQ": query sequence length
#  (5) "XT": number of continues reads, always 1 for blasr


## human readable aln fmt
# blasr winner.fa genome.fa  -m 0 

## XML fmt
# blasr winner.fa genome.fa  -m 2

## vulgar fmt --- represents sequence as symbols for match, mismatch, gap, insertion, del, etc -- also has info in beginning
# blasr winner.fa genome.fa  -m 3

## 13 field space delim -- qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
# blasr winner.fa genome.fa  -m 4 

## 19 field space delim -- qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
# blasr winner.fa genome.fa  -m 5



##Some code from Nick Loman for extracting info from BLASR
##
##import pysam
##import sys
## 
##samfile = pysam.Samfile(sys.argv[1], "rb")
## 
##fields = ['Name', 'QueryLen', 'AlignLen', 'NumMismatches']
##print "\t".join(fields)
## 
##for read in samfile:
##        t = dict(read.tags)
## 
##        results = []
##        results.append(read.qname.ljust(25))
##        results.append(t['XQ'])
##        results.append(t['XL'])
##        results.append(t['NM'])
## 
##        print "\t".join([str(r) for r in results])

