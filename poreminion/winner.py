from poretools import *
import sys

## Nov 4, 2014
## TODO for pipeline Id want both the "details" and the "fa" files
## This would require a --saveas option -- and it will save both.
## Might also want it to save the fastas to their own files when --each

#logging
import logging
logger = logging.getLogger('poreminion')


def run(parser, args):
	longest_size = 0
	longest_size_2d = 0
	longest_size_template = 0
	longest_size_complement = 0
	longest_read = None
	longest_read_2d = None
	longest_read_template = None
	longest_read_complement = None
	if args.type == 'each':
                each=True
                args.type='all'
        else:
                each=False
	for fast5 in Fast5FileSet(args.files):
		fas = fast5.get_fastas(args.type)

		for fa in fas:
                        readtype = fa.name.split()[0].split("_")[-1]
                        readlen = len(fa.seq)
                        if each:
                                if readtype == "template" and readlen > longest_size_template:
                                        longest_size_template = readlen
                                        longest_read_template = fa
                                elif (readtype == "2D" or readtype == "twodirections") and readlen > longest_size_2d:
                                        longest_size_2d = readlen
                                        longest_read_2d = fa
                                elif readtype == "complement" and readlen > longest_size_complement:
                                        longest_size_complement = readlen
                                        longest_read_complement = fa
                        else:
                                if fa and len(fa.seq) > longest_size:
                                        longest_size = len(fa.seq)
                                        longest_read = fa

		fast5.close()
##	logger.info("Wow, it's a whopper: your longest read is %d bases." % (longest_size,))
	if args.details:
                if each:
                        if longest_read_2d:
                                print ("\t").join([longest_read_2d.name.split()[0], str(longest_size_2d)])
                        if longest_read_template:
                                print ("\t").join([longest_read_template.name.split()[0], str(longest_size_template)])
                        if longest_read_complement:
                                print ("\t").join([longest_read_complement.name.split()[0], str(longest_size_complement)])
                else:
                        print ("\t").join([longest_read.name.split()[0], str(longest_size)])
        else:
                if each:
                        if longest_read_2d:
                                print longest_read_2d
                        if longest_read_template:
                                print longest_read_template
                        if longest_read_complement:
                                print longest_read_complement
                else:
                        print longest_read

