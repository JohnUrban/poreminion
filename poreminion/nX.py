#nX.py
import logging
logger = logging.getLogger('poreminion')
logger.setLevel(logging.INFO)

def NX(l, x=[25,50,75]):
        """
        Returns NX for all x for a list of numbers l.
        Default: N25, N50, N75
        Assumes all values in list x are between 0 and 100.
        Interpretation: When NX = NX_value, X% of data (in bp) is contained in reads at least NX_value bp long.
        """
	if isinstance(l, list) and isinstance(x, list):
		l = sorted(l)
		x = sorted(x)
		total = sum(l)
                nxsum = 0
                nxvalues = {e:0 for e in x}
		for e in x:
                        xpct = total*e/100.0
                        while nxsum < xpct and l:
                                nxsum += l[-1]
                                lastsize = l.pop()
                        nxvalues[e] = lastsize
                return nxvalues

	else:
		return None
	    
def run(parser, args):

    if args.inputfile:
        assert args.colnum > 0
        l = []
        try:
            for line in open(args.inputfile):
                line = line.strip().split()
                if line: ## buffer against empty lines (e.g. at end of file)
                    l.append(float(line[args.colnum-1]))
        except IndexError:
            logger.error("This file does not have this many columns. Please provide correct tab-delimited file or specify correct column.")
            quit()
    elif args.cmdline:
        l = [float(e) for e in args.cmdline.split(",")]
        
##    if args.pctreadsgtx:
##            print pctReadsGtX(l, float(args.pctreadsgtx))
##    elif args.pctdatagtx:
##            print pctDataGtX(l, float(args.pctdatagtx))
##    else: ## NX operation default
        
    l.sort()
    x = [float(e) for e in args.x.split(",")]
    nxvalues = NX(l,x)
    for e in x:
        print "N%s\t%d" % (str(e), nxvalues[e])
