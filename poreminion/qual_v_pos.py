from poretools.Fast5File import *
from collections import defaultdict
import numpy as np
import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
## matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

#logging
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)


def plot_out(args, plottype):
    '''plottype is string to be added to base filename '''
    if args.saveas is not None:
        logger.info("Writing plot to file...")
        plot_file = args.saveas.split(".")
        plot_file.insert(-1, plottype)
        plot_file = (".").join(plot_file)
        if plot_file.endswith(".pdf") or plot_file.endswith(".jpg"):
                plt.savefig(plot_file)
        else:
                logger.error("Unrecognized extension for %s! Try .pdf or .jpg" % (plot_file))
                sys.exit()

    else:
            logger.info("Showing plot...")
            plt.show()

def run(parser, args):
    """ returns dictionary with keys=positions, values=lists of qual scores for that position"""
    qualpos = defaultdict(list)
    zpos = defaultdict(list)
    bin_width = args.bin_width
    read_lengths = list()
    mean_scores = list()
##    if args.qualVsLen:
##        x = list()
##        y = list()
    for fast5 in Fast5FileSet(args.files):
        if args.start_time or args.end_time:
                read_start_time = fast5.get_start_time()
                read_end_time = fast5.get_end_time()
                if args.start_time and args.start_time > read_start_time:
                        fast5.close()
                        continue
                if args.end_time and args.end_time < read_end_time:
                        fast5.close()
                        continue

        fqs = fast5.get_fastqs(args.type)
        if args.high_quality:
                if fast5.get_complement_events_count() <= \
                   fast5.get_template_events_count():
                        fast5.close()
                        continue
        
        for fq in fqs:
                if fq is None or len(fq.seq) < args.min_length or len(fq.seq) > args.max_length:			
                        continue
                ctr = 0
##                if args.zscore:
                qscores = []
                for q in fq.qual:
                    qscores.append(ord(q)-33)
                qmean = np.mean(qscores)
                qSD = np.std(qscores)
                qLen = len(qscores)
                for q in qscores:
                    ctr += 1
                    zpos[1+int(ctr//bin_width)].append((q-qmean)/qSD)
                    qualpos[1+int(ctr//bin_width)].append(q)
                read_lengths.append(qLen)
                mean_scores.append(qmean)
##                elif args.qualVsLen:
##                    qscores = []
##                    for q in fq.qual:
##                        qscores.append(ord(q)-33)
##                    qmean = np.mean(qscores)
##                    qLen = len(qscores)
##                    x.append(qLen)
##                    y.append(qmean)
##                    
##                else:
##                    for q in fq.qual:
##                        ctr += 1
##                        qualpos[1+int(ctr//bin_width)].append(ord(q)-33)
##
        fast5.close()

## qual v readlen scatter
    if read_lengths:
        logger.info("Constructing quality vs. readlen scatter plot...")
        plt.scatter(read_lengths, mean_scores)
        plt.xlabel("Read Length (bp)")
        plt.ylabel("Mean Quality score")
        plot_out(args, "meanqual-Vs-readlen")
        plt.close()
    else: logger.warn("No reads that meet criteria: cannot construct quality vs. readlen scatter plot...")

#qual vs. pos boxplot
    if qualpos.keys():
        logger.info("Processing data...")
        data = [qualpos[e] for e in sorted(qualpos.keys())]
        logger.info("Constructing quality vs. position box plot...")
        plt.boxplot(data)
        xdetail = " (" + str(bin_width) + " bp bins)"
        plt.xlabel("Bin number in read" + xdetail)
        plt.ylabel("Quality score")
        plt.xticks(rotation=65, fontsize=8)
        plot_out(args, "qual-Vs-pos")
        plt.close()
    else: logger.warn("No reads that meet criteria: cannot construct quality vs. position box plot...")


#zscore vs. pos boxplot
    if zpos.keys():
        logger.info("Processing data...")
        data = [zpos[e] for e in sorted(zpos.keys())]
        logger.info("Constructing quality-Z-score vs. position box plot...")
        plt.boxplot(data)
        xdetail = " (" + str(bin_width) + " bp bins)"
        plt.xlabel("Bin number in read" + xdetail)
        plt.ylabel("Quality Z-score")
        plt.xticks(rotation=65, fontsize=8)
        plot_out(args, "qualZscore-Vs-pos")
        plt.close()
    else: logger.warn("No reads that meet criteria: cannot construct quality-Zscore  vs. position scatter plot...")



