import poretools.Fast5File as Fast5File
from collections import defaultdict
import pandas
import matplotlib.pyplot as plt

##TODO -- allow it to select all, just temp, just comp

def run(parser, args):
    """ returns dictionary with keys=positions, values=lists of qual scores for that position"""
    qualpos = defaultdict(list)
    bin_width = args.bin_width
    
    for fast5 in Fast5File.Fast5FileSet(args.files):
##            fq = fast5.get_fastq()
##            if fq is not None:
##                ctr = 0 
##                for q in fq.qual:
##                    ctr += 1
##                    qualpos[1+int(ctr//bin_width)].append(ord(q)-33)
##            fast5.close()
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
                for q in fq.qual:
                    ctr += 1
                    qualpos[1+int(ctr//bin_width)].append(ord(q)-33)

        fast5.close()

##    for q in sorted(qualpos.keys()):
##            print '\t'.join(str(e) for e in [q, qualpos[q]])
    print 'making data'
    data = [qualpos[e] for e in sorted(qualpos.keys())]
##    data = data[:10]
    print 'boxing'
    plt.boxplot(data)
    plt.xlabel("Position in read")
    plt.ylabel("Quality score")
    plt.xticks(rotation=65, fontsize=8)
    print 'showing'
##    plt.show()

    plt.savefig("qualpos.pdf")


