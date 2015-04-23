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
    if args.command == 'uncalled':
        import findUncalled as submodule
    elif args.command == 'timetest':
        import findTimeErrors as submodule
    elif args.command == 'fragstats':
        import fragstats as submodule
    elif args.command == 'fragsummary':
        import fragsummary as submodule
    elif args.command == 'fragrobust':
        import robust as submodule
    elif args.command == 'nx':
        import nX as submodule
    elif args.command == 'pct2d':
        import pct2D as submodule
    elif args.command == 'has2d':
        import has2D as submodule
    elif args.command == 'numevents':
        import numevents as submodule
    elif args.command == 'events':
        import get_events as submodule
    elif args.command == 'info':
        import info as submodule
    elif args.command == 'dataconc':
        import dataconc as submodule
    elif args.command == 'qualpos':
        import qual_v_pos as submodule
    elif args.command == 'kmer':
        import kmer as submodule
    elif args.command == 'kmerplot':
        import kmerplot as submodule
    elif args.command == 'kmerdiff':
        import kmerdiff as submodule
##    elif args.command == 'align':
##        import align as submodule
    elif args.command == 'winner':
        import winner as submodule
    elif args.command == 'qualdist':
        import qualdist as submodule


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
    parser = argparse.ArgumentParser(prog='poreminion', formatter_class=argparse.RawTextHelpFormatter)#ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed poreminion version",
                        action="version",
                        version="%(prog)s " + str(poreminion.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################



    ##########
    # find uncalled (not basecalled) files
    ##########
    parser_uncalled = subparsers.add_parser('uncalled',
                                        help='Find Fast5 files that were not base-called.')
    parser_uncalled.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_uncalled.add_argument('--outprefix', "-o",
                               type=str, required=True,
                              help='Uses this as basename for the following output files: (1) list of files not basecalled because template events not found, (2) list of files not basecalled because too few events found, (3) list of files not basecalled because too many events found. (4) event stats on each.')
    parser_uncalled.add_argument('--move', "-m",
                               action='store_true', default=False,
                              help='''If specified, will move each non-basecalled file type to an approp labeled dir
                                        inside same dir that has the dir reads with reads in it (e.g. downloads --> pass,
                                        downloads --> fail, downloads --> "notemplate", etc).
                                        Still writes out stats file.''')
    parser_uncalled.set_defaults(func=run_subtool)


    ##########
    # findTimeErrors
    ##########
    parser_timetest = subparsers.add_parser('timetest',
                                        help='Find Fast5 files that have event times that are earlier than event times before it suggesting malfunction/erroneous read.')
    parser_timetest.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_timetest.add_argument('--outprefix', "-o",
                               type=str, default=False,
                              help='Uses this as basename for file containing list of files with time errors.')
    parser_timetest.add_argument('--move', "-m",
                               action='store_true', default=False,
                              help='''If specified, will move files with time error dir labeled time_errors
                                        inside same dir that has the dir with reads in it (e.g. downloads --> pass,
                                        downloads --> fail, downloads --> "time_errors", etc).
                                        Still writes out list file above.''')
    parser_timetest.add_argument('--verbose', "-v",
                               action='store_true', default=False,
                              help='''Will print to stderr info about how far along it is in process.''')
    parser_timetest.set_defaults(func=run_subtool)




    ##########
    # fragstats
    ##########
    parser_fragstats = subparsers.add_parser('fragstats',
                                        help='''Run this on set of base-called fast5 files.
Returns tab-delimited table with columns:
1 = readname,
2 = estimated molecule/fragment size,
3 = number input events,
4 = if complement detected,
5 = if 2D detected,
6 = num template events,
7 = num complement events,
8 = length of 2D sequence,
9 = length of template sequence,
10 = length of complement sequence,
11 = mean qscore of 2D sequence,
12 = mean qscore of template sequence,
13 = mean qscore of complement,
14 = ratio of number template events to number complement events,
15 = channel number molecule traversed
16 = heat sink temperature while molecule traversed
17 = num called template events (after events pruned during base-calling)
18 = num called complement events (after events pruned during base-calling)
19 = num skips in template
20 = num skips in complement
21 = num stays in template
22 = num stays in complement
23 = strand score template
24 = strand score complement
25 = num stutters in template
26 = num stutters in complement

If --extensive used:
27 = starttime,
28 = endtime,
29 = slope across all events,
30 = mean duration across all events,
31 = median duration across all events,
32 = sd of all event durations,
33 = min event duration,
34 = max event duration,
35-40 = num temp events with 0,1,2,3,4,5 moves from base-caller,
41-46 = num comp events with 0,1,2,3,4,5 moves from base caller.

If --checktime used:
Final column = 0 or 1 for no/yes there is a time error present.

Estimates molecule/fragment size in the following way.
If has 2D, molecule size is the length of 2D read.
If template only, molecule size is the length of template read.
If template and complement, but no 2D, molecule size is length of the longer read between template and complement.
Molecule size allows calculation of total non-redundant data.
This is the sum of unique molecule lengths rather than summing all read types from each molecule.
From the molecule sizes, the "Molecule N50" can be computed using the nx subcommand on the fragstats file and specifying colum 2.
                                                    ''')
    parser_fragstats.add_argument('--extensive', "-e",
                               action="store_true", required=False, default=False,
                              help='''This tacks a number of fields on at the end of the regular frag stats that requires much much more computation time.
                                    The additional fields are: 16=starttime of all events, 17=endtime of all events, 18=slope of all events,
                                    19=mean duration across all events, 20=median duration across all events, 21=sd of all event durations, 22=min event duration, 23=max event duration,
                                    24-29=num temp events with 0,1,2,3,4,5 moves from base-caller, 20-35=num comp events with 0,1,2,3,4,5 moves from base caller.
                                    ''')
    parser_fragstats.add_argument('--checktime', "-t",
                               action="store_true", required=False, default=False,
                              help='''This tacks on timetest info (search for time errors in start times) as the last field
                                    --> 0 or 1 for no/yes there is a time error present. Adds considerable computation time.
                                    If used with --extensive, will take even more time than that alone.''')

    parser_fragstats.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    
    parser_fragstats.set_defaults(func=run_subtool)


    ##########
    # fragsummary
    ##########
    parser_fragsummary = subparsers.add_parser('fragsummary',
                                        help='''To summarize fragstats, use this with a tab-delimited, fragstats table file (output of fragstats subcommand).''')
    parser_fragsummary.add_argument('--fragfile', "-f",
                               type=str, required=True,
                              help='''Specify path to the fragstats table file (output of fragstats subcommand).
                                    ''')    
    parser_fragsummary.add_argument('--extensive', "-e",
                               action="store_true", required=False, default=False,
                              help='''Use this flag if the fragstats file was generated with -e/--extensive option.
                                    ''')
    parser_fragsummary.add_argument('--checktime', "-t",
                               action="store_true", required=False, default=False,
                              help='''Use this flag if the fragstats file was generated with -t/--checktime option.''')

    
    parser_fragsummary.set_defaults(func=run_subtool)


    ##########
    # fragsort/plot
    ##########


    

    ##########
    # nX
    ##########
    parser_nx = subparsers.add_parser('nx',
                                        help='Computes N50 or NX values on columns of a file or from comma-separated list.')

    parser_nx_input = parser_nx.add_mutually_exclusive_group(required=True)
    parser_nx_input.add_argument('-i', "--inputfile",
                       type= str, default=False,
                       help='''Input file.''')
    parser_nx_input.add_argument('--cmdline', '-c',
                       type= str, default=False,
                       help='''Input list of numbers on cmd line (comma-separated) -- e.g. -c 3,5,10,30,11 ''')

    parser_nx.add_argument('-k', "--colnum",
                       type=int, default=1,
                       help='''Column number (1-based) to compute n50 on from Input file. Default is first column.''')

    parser_nx.add_argument('-x', "--x",
                       type=str, default="25,50,75",
                       help='''Give comma-separated X values for NX function -- i.e. 50 for N50. Default=25,50,75''')

##    parser_nx.add_argument('-pctdatagtx',
##                       type=str, default=False,
##                       help='''Instead of NX values, return pct of data from lengths greater than X. Provide X with this flag.''')
##    parser_nx.add_argument('-pctreadsgtx',
##                       type=str, default=False,
##                       help='''Instead of NX values, return pct of items (reads, contigs, etc) in list greater than X. Provide X with this flag.''')

    parser_nx.set_defaults(func=run_subtool)   


    ##########
    # robust
    ##########
    parser_robust = subparsers.add_parser('fragrobust',
                                        help='''Looks at fragsizes in fragstats. Sees what percent of fragsizes are "robust" to all sequence lengths from same molecule.''')
    parser_robust.add_argument('--fragfile', "-f",
                               type=str, required=False, default=False,
                              help='''Specify path to the fragstats table file (output of fragstats subcommand).
                                    ''')    
    parser_robust.add_argument('--message', "-m",
                               action="store_true", required=False, default=False,
                              help='''Use this flag print caution message on this metric. If used with -f, it is header message.
                                    ''')
    parser_robust.set_defaults(func=run_subtool)

    ##########
    # pct 2D
    ##########
    parser_pct2d = subparsers.add_parser('pct2d',
                                        help='Get the proportion of reads that have a 2D read')
    parser_pct2d.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_pct2d.set_defaults(func=run_subtool)


    ##########
    # has 2D
    ##########
    parser_has2d = subparsers.add_parser('has2d',
                                        help='Prints 2 columns: filename, has2D =  True/False')
    parser_has2d.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_has2d_filter = parser_has2d.add_mutually_exclusive_group()
    parser_has2d_filter.add_argument('--only2d', "-2",
                               action='store_true', default=False,
                              help='''If specified, will only print out files that have 2D -- no True/False column.''')
    parser_has2d_filter.add_argument('--no2d', "-0",
                               action='store_true', default=False,
                              help='''If specified, will only print out files that do not have 2D -- no True/False column.''')
 
    parser_has2d.set_defaults(func=run_subtool)





    ##########
    # get num events
    ##########
    parser_numevents = subparsers.add_parser('numevents',
                                        help='Print 2 column list of file and number of input events in file.')
    parser_numevents.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_numevents.set_defaults(func=run_subtool)


    ##########
    # get_events
    ##########

    parser_get_events = subparsers.add_parser('events',
                                        help='''Look at events inside raw and basecalled fast5 files. ''')
    parser_get_events.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    
    parser_get_events_filetype = parser_get_events.add_mutually_exclusive_group(required=False)
    parser_get_events_filetype.add_argument('-r', '--raw', action='store_true', default=False)
    parser_get_events_filetype.add_argument('-b', '--basecalled', action='store_true', default=False)

    parser_get_events.add_argument("-t", "--type", choices=['input', 'template', 'complement'], default="input",
                               help='''What events should be returned? Specify: input, template, complement. Default: input.
                                    Template and complement events can only be specified from basecalled fast5 files.''')

    parser_get_events.add_argument("-H", "--header", action="store_true", default=False,
                               help='''Adds header line to top with a "#" at beginning of line.''')
    
    parser_get_events.set_defaults(func=run_subtool)



    ##########
    # info
    ##########

    parser_info = subparsers.add_parser('info',
                                        help='''Get info about run and, if file is basecalled, basecalling. ''')
    parser_info.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    

    parser_info.add_argument("-b", '--basic', action="store_true", default=False, help='''Some basic info.''')
    parser_info.add_argument("-a", '--all', action="store_true", default=False, help='''All info.''') 
    parser_info.set_defaults(func=run_subtool)



    ##########
    # data_conc (data concentration plot)
    ##########
    parser_dataconc = subparsers.add_parser('dataconc',
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
    parser_qualpos.add_argument('--zscore',
                              default=False,
                              action='store_true',
                              help='For each read, normalize each bucket score to the mean and stdDev of all scores in read. Z = (bucketScore-mean)/stdDev')
    parser_qualpos.add_argument('--qualVsLen',
                          default=False,
                          action='store_true',
                          help='Scatter plot mean score (y-axis) vs. read length (x-axis)')
    parser_qualpos.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save the plot to a file named filename.extension (e.g. pdf, jpg)''',
                             default=None)

    parser_qualpos.set_defaults(func=run_subtool)




    ##########
    # qualdist
    ##########
    parser_qualdist = subparsers.add_parser('qualdist',
                                        help='''Get the qual score composition of a set of FAST5 files.
This tool is from poretools, but poreminion allows you to select the type of read.''')
    parser_qualdist.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_qualdist.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev'],
                              default='all',
                              help='Which type of reads should be analyzed? Def.=all, choices=[all, fwd, rev, 2D, fwd,rev]. Is mutually exclusive with --one-read-per-molecule.')
 
    parser_qualdist.set_defaults(func=run_subtool)





    ##########
    # kmerCounting
    ##########
    parser_kmer = subparsers.add_parser('kmer',
                                        help='Count kmers in reads or reference.')
    parser_kmer.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_kmer.add_argument('-k', '--kmersize',
                              dest='k',
                              default=5,
                              type=int,
                              help=('Kmer size. Default = 5. Sizes 1-7 work well with kmerplot on regular Mac OS. Up to 10 is possible. After that it might require too much memory for kmerplot on regular Mac OS.'))
    parser_kmer.add_argument('--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file instead of fast5dir/.
                                    While min and max length arguments remain meaningful for fasta files, the following arguments do not: start time, end time, high quality, type, single read per molecule.'''))
    parser_kmer.add_argument('--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file instead of fast5dir/.
                                    While min and max length arguments remain meaningful for fastq files, the following arguments do not: start time, end time, high quality, type, single read per molecule.'''))
    parser_kmer.add_argument('--rev-comp',
                              dest='rev_comp',
                              default=False,
                              action="store_true",
                              help='''Created to be used with --fasta and --fastq options.
                                    When creating kmer counts, it counts both the fwd and reverse complement kmers.
                                    For now, it does nothing when used with fast5 dirs (minION data files).''')

##    parser_kmer_output = parser_kmer.add_mutually_exclusive_group()
##    parser_kmer_output.add_argument('-t', '--table',
##                              dest='table',
##                              default=True,
##                              action='store_true',
##                              help=('''Output option: report tab-delimited table of kmer, count, and proportion of all kmers seen.
##                                    Default = True (to stdout). Use --saveas to specify file to save to.'''))
##    parser_kmer_output.add_argument('-p', '--plot',
##                              dest='plot',
##                              default=False,
##                              action='store_true',
##                              help=('''Output option: show or write out plot.
##                                    Default = False (to stdout). Use --saveas to specify file to save to.'''))
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

    
    ##########
    # kmerplotting
    ##########
    parser_kmerplot = subparsers.add_parser('kmerplot',
                                        help='Plot kmer counts in reads or reference.')
##    parser_kmerplot.add_argument('files', metavar='FILES', nargs='+',
##                             help='The input FAST5 files.')

    parser_kmerplot.add_argument('-t1', '--kmer-count-in-reads',
                             dest='table1',
                             type=str,
                             help='''Provide path to file with kmer count table from reads (or any kmer count table).
                                    This argument is required and when used alone, just generates a bar plot of kmer counts.''',
                             default=None)
    
    parser_kmerplot.add_argument('-t2', '--kmer-count-in-reference',
                             dest='table2',
                             type=str,
                             help='''Provide path to file with kmer count table from reference sequence (or any second kmer count table).
                                    This argument is not required and if used, results in a scatterplot of the 2 kmer count tables.''',
                             default=None)
    parser_kmerplot.add_argument('--matplotlib',
                             dest='mpl',
                             action='store_true',
                             help='''Temp option: plot in matplotlib''',
                             default=False)
    parser_kmerplot.add_argument('--ggplot2',
                             dest='gg',
                             action='store_true',
                             help='''Temp option: plot in ggplot2''',
                             default=False)
    parser_kmerplot.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" where extension can be only pdf and jpg for now.''',
                             default=None)
    parser_kmerplot.set_defaults(func=run_subtool)
    

    ##########
    # kmer diff abundance
    ##########
    parser_kmerdiff = subparsers.add_parser('kmerdiff',
                                        help='Get fold-enrichment values of kmers in reads vs reference.')

    parser_kmerdiff.add_argument('-t1', '--kmer-count-in-reads',
                             dest='table1',
                             type=str,
                             help='''Provide path to file with kmer count table from reads (or any kmer count table).
                                    This argument is required and when used alone, just generates a bar plot of kmer counts.''',
                             default=None)
    
    parser_kmerdiff.add_argument('-t2', '--kmer-count-in-reference',
                             dest='table2',
                             type=str,
                             help='''Provide path to file with kmer count table from reference sequence (or any second kmer count table).
                                    This argument is not required and if used, results in a scatterplot of the 2 kmer count tables.''',
                             default=None)
    parser_kmerdiff.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" where extension can be only pdf and jpg for now.''',
                             default=None)
    parser_kmerdiff.add_argument('-bcv', '--square-root-dispersion',
                             dest='bcv',
                             type=float,
                             help='''When there are no replicates in edgeR, dispersion must be determined by the user.
                                    The default is 0.2. Other values to try could be 0.01-0.4 (or any).
                                    p-values will be sensitive to choice of bcv. Fold change will not.''',
                             default=0.2)

    parser_kmerdiff.add_argument('--volcano',
                             dest='volcano',
                             type=str,
                             help='''If you want the analysis to generate a volcano plot,
                                    (log(fold change) vs. -log10(pvalue)), then use this flag
                                    and provide the name and extension of volcano plot file (e.g. volcano.jpg).''',
                             default=None)
    parser_kmerdiff.add_argument('--smear',
                             dest='smear',
                             type=str,
                             help='''If you want the analysis to generate a smear plot,
                                    (log(fold change) vs. log(CPM)), then use this flag
                                    and provide the name and extension of smear plot file (e.g. smear.jpg).''',
                             default=None)
    parser_kmerdiff.add_argument('--nt-content',
                             dest='nt_content',
                             type=str,
                             help='''If you want the analysis to generate a table analyzing the nucleotide content
                                    of kmers >= abs(fold change) and pval <= p, then use this flag with those values as in these
                                    examples: (a) --nt-content fc:2.5,p:0.001 (b) --nt-content fc:2,fdr:0.1)''',
                             default=None)
    
    parser_kmerdiff.set_defaults(func=run_subtool)



    ##########
    # winner -- adds "each" and "details" functionalities to poretools winner
    ##########
    parser_winner = subparsers.add_parser('winner',
                                        help='''Get the longest read from a set of FAST5 files.
Similar to poretools winner, only allows type=each and offers a details only option.''')
    parser_winner.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_winner.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev', 'each'],
                              default='all',
                              help='''Which type of FASTA entries should be reported? Def.=all.
                                    Choices: 'all', 'fwd', 'rev', '2D', 'fwd,rev', 'each'.
                                    'each' will give longest for each 2D, fwd, rev.''')
    parser_winner.add_argument('--details',
                               action="store_true",
                              default=False,
                              help='If set, it will only print details: readname, length')
    parser_winner.set_defaults(func=run_subtool)



    ##########
##    # viz time
##    ##########
##    parser_numevents = subparsers.add_parser('viztime',
##                                        help='Visualize the relative start time across events - e.g. to visualize a time error or lack thereof.')
##    parser_numevents.add_argument('files', metavar='FILES', nargs='+',
##                               help='The input FAST5 files.')
##    parser_numevents.set_defaults(func=run_subtool)

    ##########
    # alignment
    ##########
##    parser_align = subparsers.add_parser('align',
##                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Performs alignments -- returns alignments, stats, plots')
##    align_subparsers = parser_align.add_subparsers(title='[align-commands]', dest='align_command', parser_class=ArgumentParserWithDefaults)
##    parser_align_fitting = align_subparsers.add_parser('fitting', help="fitting of 1 DNA seq to another")
##    parser_align_blasr = align_subparsers.add_parser('blasr', help="BLASR")
##    parser_align_blastn = align_subparsers.add_parser('blastn', help="BLASTN")
##    parser_align_last = align_subparsers.add_parser('last', help="LAST")
##    ## BLASR
##    parser_align_blasr.add_argument("--plot")
##    parser_align_blasr.add_argument("--blasr_file", type=str,default=None)
##    parser_align.add_argument('--sequence',
##                             dest='sequence',
##                             type=str,
##                             help='''Provide a sequence < min read length. Default is the KanRR fragment: CGTACGCTGCAGGTCG''',
##                             default='CGTACGCTGCAGGTCG') ##shortest version of KanR fragment that maximized scores on some pilot reads
##    parser_align.add_argument('-m', '--multiple-sequences',
##                             dest='multiple_sequences',
##                             type=str, default=None,
##                             help='''Provide a path to a file with 1 sequence per line.
##                                    For each read in the fastx file, it will report the fitting alignment for the
##                                    sequence in this file with the best fitting aln.
##                                    Each time it encounters a score that ties the current max score, it exchanges the older fiting aln
##                                    info for the new fitting aln info with a 50%% probability.
##                                    This way there is a random assignment of the best barcode.
##                                    Use --all-scores instead to get an output with all max scores and barcodes returned.''')
##    parser_align.add_argument('-w', '--with-read-names',
##                             dest='with_read_names', action="store_true",
##                             default=False,
##                             help='''If set, will print "readname, startPosInRead, fitAlnScore, fitAlnScore/queryLen";
##                                else just "startPosInRead,fitAlnScore, fitAlnScore/queryLen".
##                                Start position is in pythonese (0-based).''')
##    parser_align.add_argument('-e', '--with-edit-distances',
##                             dest='with_edit_distances', action="store_true",
##                             default=False,
##                             help='''If set, edit dist will be incl in output''')
##    parser_align.add_argument('-a', '--with-aln-seqs',
##                             dest='with_aln_seqs', action="store_true",
##                             default=False,
##                             help='''If set, the aligned versions of sequences 1 (read) and 2 (provided) will be printed.''')
##    parser_align.add_argument('-r', '--random-sequence',
##                             dest='random_sequence', type=int,
##                             default=False,
##                             help='''Provide integer for random sequence length. This option overrides --sequence.''')
##    parser_align_seqtransform = parser_align.add_mutually_exclusive_group()
##    parser_align_seqtransform.add_argument("-c", "--complement", action="store_true", default=False,
##                                                 help=''' Use complement of provided sequence -- right now only works on single seq.
##                                                            e.g. AACC -> TTGG''')
##    parser_align_seqtransform.add_argument("-rc", "--reverse_complement", action="store_true", default=False,
##                                                 help=''' Use reverse complement of provided sequence -- right now only works on single seq.
##                                                            e.g. AACC -> GGTT''')
##    parser_align_seqtransform.add_argument("-rs", "--reverse_sequence", action="store_true", default=False,
##                                                 help=''' Use reverse sequence of provided sequence -- right now only works on single seq.
##                                                            e.g. AACC -> CCAA''')
##
##    parser_align.set_defaults(func=run_subtool)
##

                                     


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
