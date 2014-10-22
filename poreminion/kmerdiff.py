## TODO: Clean up this file. Do away with unused fxns. Do away with need for "Total". etc.
##      Break run edger up -- have it return the DF, have each operation (e.g. volcano plot) be its own fxn
##      Have run() be the only fxn that handles parser/arg

"""Calculate differentially represented kmers using EdgeR from bioconductor.

http://bioconductor.org/packages/2.5/bioc/html/edgeR.html

This analysis is derived from Brad Chapman's example at:
https://github.com/chapmanb/bcbb/blob/master/stats/count_diffexp.py
"""
import os
import sys
import csv
from collections import defaultdict
import numpy
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from kmerplot import ensureEqualKmerSets, kmerDictToPlotData, readInTwoKmerTables, totalKmerCount
import matplotlib.pyplot as plt
from math import log10, log
#logging
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)


def log2(x):
    return log(x)/float(log(2))

def basecounter(kmer, basedict):
    '''basedict is a defaultdict(int) -- provide empty or non-empty one for update'''
    for b in kmer.upper():
        basedict[b] += 1
    return basedict

def run_edger(data, groups, sizes, genes, args):
    """Call edgeR in R and organize the resulting differential expressed genes.
    """
    robjects.r('''
        library(edgeR, quietly = TRUE)
    ''')

    robjects.r('''
            f <- function(data, groups, genes, bcv){
            y <- DGEList(counts=data, group=groups, genes=as.matrix(genes, ncol=1))
            y <- calcNormFactors(y)
            y <- calcNormFactors(y)
            et <- exactTest(y, dispersion=bcv^2)
            return(et)
            }
            ''')
    r_f = robjects.globalenv['f']
    dgelist = r_f(data, groups, genes, args.bcv)
    tag_table = robjects.r.topTags(dgelist, n=len(genes))[0]
    p = []
    fc = []
    k = []
    cpm = []
    fdr = []
    if args.saveas is not None:
        logger.info("columns: kmerNumber, kmer, logFC, logCPM, Pvalue, FDR")
        fout = open(args.saveas, 'w')
        for e in tag_table.iter_row():
            out = ("\t").join(str(e).split("\n")[1].split())
            var = out.split("\t")
            p.append(float(var[4]))
            fc.append(float(var[2]))
            k.append(var[1])
            cpm.append(float(var[3]))
            fdr.append(float(var[5]))
            fout.write(out+"\n")
    else:
##        print tag_table
        for e in tag_table.iter_row():
            out = ("\t").join(str(e).split("\n")[1].split())
            var = out.split("\t")
            p.append(float(var[4]))
            fc.append(float(var[2]))
            k.append(var[1])
            cpm.append(float(var[3]))
            fdr.append(float(var[5]))
            print out
    if args.volcano:
        logger.info("Generating Volcano plot.")
        x = [e for e in fc]
        y = [-1*log10(e) for e in p]
        plt.scatter(x=x, y=y, s=5)
        for i in range(len(k)):
            plt.annotate(k[i], (x[i], y[i]), size='xx-small')
        plt.xlabel("log2(Fold Change)")
        plt.ylabel("-log10(p-value)")
##        plt.show()
        plot_file = args.volcano
        if plot_file.endswith(".pdf") or plot_file.endswith(".jpg"):
                plt.savefig(plot_file)
        plt.close()
    if args.smear:
        logger.info("Generating smear plot.")
        x = [e for e in fc]
        y = [e for e in cpm]
        plt.scatter(x=x, y=y, s=5)
        for i in range(len(k)):
            plt.annotate(k[i], (x[i], y[i]), size='xx-small')
        plt.xlabel("log2(Fold Change)")
        plt.ylabel("log(2CPM)")
##        plt.show()
        plot_file = args.smear
        if plot_file.endswith(".pdf") or plot_file.endswith(".jpg"):
                plt.savefig(plot_file)
        plt.close()
    if args.nt_content:
        proceed = True
        filt = {'p':float('inf'), 'fdr':float('inf'), 'fc':1e-100}
        params = args.nt_content.split(",")
        for e in params:
            x,y = e.split(":")
            y = float(y)
            try:
                filt[x] = y
            except KeyError:
                logger.error("KeyError: provide only 'p', 'fc', and/or 'fdr' in a comma-separated k:v,k:v format")
                proceed = False
                break
        if proceed:
            basedictFE = defaultdict(int) #enriched
            basedictFD = defaultdict(int) #depleted
            for i in range(len(k)):
##                print abs(fc[i]) >= filt['fc'], p[i] <= filt['p'], fdr[i] <= filt['fdr']
                if abs(fc[i]) >= log2(filt['fc']) and p[i] <= filt['p'] and fdr[i] <= filt['fdr']:
                    if fc[i] >= log2(filt['fc']):
                        basedictFE = basecounter(k[i], basedictFE)
                    else:
                        basedictFD = basecounter(k[i], basedictFD)
            totalFE = float(sum(basedictFE.values()))
            totalFD = float(sum(basedictFD.values()))
            print "Enriched Kmers:"
            for k,v in sorted(basedictFE.iteritems()):
                print ("\t").join([k,str(v), str(v/totalFE)])
            basedictFE['AT'] = basedictFE['A'] + basedictFE['T']
            basedictFE['GC'] = basedictFE['G'] + basedictFE['C']
            print ("\t").join(['AT',str(basedictFE['AT']), str(basedictFE['AT']/totalFE)])
            print ("\t").join(['GC',str(basedictFE['GC']), str(basedictFE['GC']/totalFE)])
            print "Depleted Kmers:"
            for k,v in sorted(basedictFD.iteritems()):
                print ("\t").join([k,str(v), str(v/totalFD)])
            basedictFD['AT'] = basedictFD['A'] + basedictFD['T']
            basedictFD['GC'] = basedictFD['G'] + basedictFD['C']
            print ("\t").join(['AT',str(basedictFD['AT']), str(basedictFD['AT']/totalFD)])
            print ("\t").join(['GC',str(basedictFD['GC']), str(basedictFD['GC']/totalFD)])
            ## TODO PWM/Logo over kmer


def get_conditions_and_genes(work_counts): 
    conditions = work_counts.keys()
    conditions.sort()
    all_genes = []
    for c in conditions:
        all_genes.extend(work_counts[c].keys())
    all_genes = list(set(all_genes))
    all_genes.sort()
    sizes = [work_counts[c]["Total"] for c in conditions]
    all_genes.remove("Total")
    return conditions, all_genes, sizes
    
def edger_matrices(work_counts):
    """Retrieve matrices for input into edgeR differential expression analysis.
    """
    conditions, all_genes, sizes = get_conditions_and_genes(work_counts) #get_conditions_and_genes
    assert len(sizes) == 2
    groups = [1, 2]
    data = []
    final_genes = []
    for g in all_genes:
        cur_row = [int(work_counts[c][g]) for c in conditions]
        if sum(cur_row) > 0:
            data.append(cur_row)
            final_genes.append(g)
    return (numpy.array(data), numpy.array(groups), numpy.array(sizes),
            conditions, final_genes)

##def read_count_file(in_file):
##    """Read count information from a simple CSV file into a dictionary.
##    """
##    counts = defaultdict(dict)
##    with open(in_file) as in_handle:
##        reader = csv.reader(in_handle)
##        header = reader.next()
##        conditions = header[1:]
##        for parts in reader:
##            region_name = parts[0]
##            region_counts = [float(x) for x in parts[1:]]
##            for ci, condition in enumerate(conditions):
##                counts[condition][region_name] = region_counts[ci]
##    return dict(counts)

##def main(count_file):
##    base, ext = os.path.splitext(count_file)
##    outfile = "%s-diffs.csv" % (base)
##    counts = read_count_file(count_file)
##    data, groups, sizes, conditions, genes = edger_matrices(counts)
##    probs = run_edger(data, groups, sizes, genes)
##    write_outfile(outfile, genes, conditions, counts, probs)

##def create_count_file(parser, args):
##    ## read in and equilibrate the 2 kmer count tables
##    kmerdict1, kmerdict2 = readInTwoKmerTables(parser, args)
##    ## get total counts
##    total1 = totalKmerCount(kmerdict1)
##    total2 = totalKmerCount(kmerdict2)
##    ## write out header
##    prefix = args.saveas.split(".")[0]
##    count_fh = "%s.count.csv" % (prefix)
##    countfile = open(count_fh, 'w')
##    ## write header
##    countfile.write("kmer,condition1,condition2\n")
##    ## write total
##    countfile.write(",".join(["Total", str(total1), str(total2)]) +"\n")
##    ## write entries
##    for e in sorted(kmerdict1.keys()):
##        out = (",").join([e, str(kmerdict1[e]), str(kmerdict2[e])])
##        countfile.write(out+"\n")
##    return count_fh
    
def make_count_dict(parser, args):
    ## read in and equilibrate the 2 kmer count tables
    kmerdict1, kmerdict2 = readInTwoKmerTables(parser, args)
    ## get total counts
    total1 = totalKmerCount(kmerdict1)
    total2 = totalKmerCount(kmerdict2)
    ## add totals to respective kmerdicts
    kmerdict1["Total"] = total1
    kmerdict2["Total"] = total2
    ## add kmdericts to defaultdict(dict) that has keys condition1 and condition2
    counts = defaultdict(dict)
    counts['condition1'] = dict(kmerdict1)
    counts['condition2'] = dict(kmerdict2)
    return dict(counts)
    
    
def run(parser, args):
##    count_file = create_count_file(parser, args)
##    base, ext = os.path.splitext(count_file)
##    outfile = "%s-diffs.csv" % (base)
##    counts = read_count_file(count_file) #dict condition:kmerdict
    counts = make_count_dict(parser, args)
    data, groups, sizes, conditions, genes = edger_matrices(counts)
    probs = run_edger(data, groups, sizes, genes, args)
##    write_outfile(outfile, genes, conditions, counts, probs)
    
