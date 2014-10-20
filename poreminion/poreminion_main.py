#!/usr/bin/env python

import os.path
import sys
import argparse
from poretools.Fast5File import *

#logger
import logging
logger = logging.getLogger('poreminion')

# poreminion imports
import poreminion.version

def run_subtool(parser, args):
    if args.command == 'data_conc':
        import dataconc as submodule
    elif args.command == 'qualpos':
        import qual_v_pos as submodule
    elif args.command == 'kmer':
        import kmer as submodule

    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
	self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def main():
    logging.basicConfig()

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='poreminion', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed poreminion version",
                        action="version",
                        version="%(prog)s " + str(poreminion.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################


    ##########
    # data_conc (data concentration plot)
    ##########
    parser_dataconc = subparsers.add_parser('data_conc',
                                        help='''Plot sum of read lengths in each bin for a given set of bins for a set of FAST5 files.
                                        This is the type of plot seen in MinKNOW while sequencing.''')
    parser_dataconc.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_dataconc.add_argument('--min-length',
                              dest='min_length',
                              default=0,
                              type=int,
                              help=('Minimum read length to be included in analysis.'))
    parser_dataconc.add_argument('--max-length',
                              dest='max_length',
                              default=1000000000,
                              type=int,
                              help=('Maximum read length to be included in analysis.'))
    parser_dataconc.add_argument('--bin-width',
                              dest='bin_width',
                              default=500,
                              type=int,
                              help=('The width of bins (default: 500 bp).'))
    parser_dataconc.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save the plot to a file named filename.extension (e.g. pdf, jpg)''',
                             default=None)
    parser_dataconc.add_argument('--cumulative',
                                action="store_true",
                             help='''For cumulative plot.''',
                             default=False)
    parser_dataconc.add_argument('--percent',
                             action="store_true",
                             help='''Plot as percentge of all data.''',
                             default=False)
    parser_dataconc.add_argument('--simulate',
                             action="store_true",
                             help='''This will randomly sample N read lengths in the size range from min to max (or according parameters set by --parameters),
                                    where N is the number of reads in the fast5 dir (or N specified with --parameters). 
                                    Then it will plot the simulation lengths. INFO about parameters used is printed so that
                                    it can be reproduced with --parameters in the future (much faster).''',
                             default=False)
    parser_dataconc.add_argument('--parameters',
                             type=str,
                             help='''--simulate by default will use N=readcount, range=min-to-max. Override this with --parameters N,min,max. e.g. --parameters 350,500,48502''',
                             default=False)
#
    parser_dataconc.add_argument('--start',
                              dest='start_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from after start timestamp')
    parser_dataconc.add_argument('--end',
                              dest='end_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from before end timestamp')
    parser_dataconc.add_argument('--high-quality',
                              dest='high_quality',
                              default=False,
                              action='store_true',
                              help='Only analyze reads with more complement events than template. Can be used with --type or --one-read-per-molecule to select a specific read type from high quality reads.')
    parser_dataconc_readfilter = parser_dataconc.add_mutually_exclusive_group()
    parser_dataconc_readfilter.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev'],
                              default='all',
                              help='Which type of reads should be analyzed? Def.=all, choices=[all, fwd, rev, 2D, fwd,rev]. Is mutually exclusive with --one-read-per-molecule.')
    parser_dataconc_readfilter.add_argument('-1', '--one-read-per-molecule',
                              dest='single_read',
                              default=False,
                              action='store_true',
                              help='''Only analyze one read per molecule in priority order: 2D -> template -> complement.
                                            That is, if there is a 2D read use that.If not, then try to use template. etc.
                                            Is mutually exclusive with --type.''')
    parser_dataconc.set_defaults(func=run_subtool)


    ##########
    # qual vs. position
    ##########
    parser_qualpos = subparsers.add_parser('qualpos',
                                        help='Get the qual score distribution over positions in reads')
    parser_qualpos.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_qualpos.set_defaults(func=run_subtool)
    parser_qualpos.add_argument('--min-length',
                              dest='min_length',
                              default=0,
                              type=int,
                              help=('Minimum read length to be included in analysis.'))
    parser_qualpos.add_argument('--max-length',
                              dest='max_length',
                              default=1000000000,
                              type=int,
                              help=('Maximum read length to be included in analysis.'))
    parser_qualpos.add_argument('--bin-width',
                              dest='bin_width',
                              default=1000,
                              type=int,
                              help=('The width of bins (default: 1000 bp).'))
    parser_qualpos.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev'],
                              default='all',
                              help='Which type of reads should be analyzed? Def.=all, choices=[all, fwd, rev, 2D, fwd,rev]')
    parser_qualpos.add_argument('--start',
                              dest='start_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from after start timestamp')
    parser_qualpos.add_argument('--end',
                              dest='end_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from before end timestamp')
    parser_qualpos.add_argument('--high-quality',
                              dest='high_quality',
                              default=False,
                              action='store_true',
                              help='Only analyze reads with more complement events than template.')
    parser_qualpos.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save the plot to a file named filename.extension (e.g. pdf, jpg)''',
                             default=None)



    ##########
    # kmerCounting
    ##########
    parser_kmer = subparsers.add_parser('kmer',
                                        help='NEW FEATURE -- NOT YET STABLE. Count kmers in reads or reference.')
    parser_kmer.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_kmer.set_defaults(func=run_subtool)
    parser_kmer.add_argument('-k', '--kmersize',
                              dest='k',
                              default=5,
                              type=int,
                              help=('Kmer size. Default = 5.'))
    parser_kmer.add_argument('--min-length',
                              dest='min_length',
                              default=0,
                              type=int,
                              help=('Minimum read length to be included in analysis.'))
    parser_kmer.add_argument('--max-length',
                              dest='max_length',
                              default=1000000000,
                              type=int,
                              help=('Maximum read length to be included in analysis.'))
    parser_kmer.add_argument('--start',
                              dest='start_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from after start timestamp')
    parser_kmer.add_argument('--end',
                              dest='end_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from before end timestamp')
    parser_kmer.add_argument('--high-quality',
                              dest='high_quality',
                              default=False,
                              action='store_true',
                              help='Only analyze reads with more complement events than template.')
    parser_kmer.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save tab-delimited kmer + counts to file.''',
                             default=None)
    parser_kmer_readfilter = parser_kmer.add_mutually_exclusive_group()
    parser_kmer_readfilter.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev'],
                              default='all',
                              help='Which type of reads should be analyzed? Def.=all, choices=[all, fwd, rev, 2D, fwd,rev]. Is mutually exclusive with --one-read-per-molecule.')
    parser_kmer_readfilter.add_argument('-1', '--one-read-per-molecule',
                              dest='single_read',
                              default=False,
                              action='store_true',
                              help='''Only analyze one read per molecule in priority order: 2D -> template -> complement.
                                            That is, if there is a 2D read use that.If not, then try to use template. etc.
                                            Is mutually exclusive with --type.''')
    parser_kmer.set_defaults(func=run_subtool)

    


    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
      args.func(parser, args)
    except IOError, e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
