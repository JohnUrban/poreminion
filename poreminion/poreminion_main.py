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
    elif args.command == 'kmerplot':
        import kmerplot as submodule
    elif args.command == 'kmerdiff':
        import kmerdiff as submodule
    elif args.command == 'events_stats':
        import events_stats as submodule
    elif args.command == 'get_events':
        import events_stats as submodule
    elif args.command == 'get_model':
        import events_stats as submodule
    elif args.command == 'get_metadata':
        import events_stats as submodule
    elif args.command == 'align':
        import align as submodule
    elif args.command == 'winner':
        import winner as submodule
    elif args.command == 'qualdist':
        import qualdist as submodule
    elif args.command == 'pct2d':
        import pct2D as submodule
    elif args.command == 'has2d':
        import has2D as submodule
    elif args.command == 'fragstats':
        import fragstats as submodule
    elif args.command == 'uncalled':
        import findUncalled as submodule
    elif args.command == 'timetest':
        import findTimeErrors as submodule
    elif args.command == 'numevents':
        import numevents as submodule

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
                                        help='Get the qual score composition of a set of FAST5 files')
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
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Count kmers in reads or reference.')
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
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. plot kmer counts in reads or reference.')
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
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. plot kmer counts in reads or reference.')

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
    # events stats
    ##########
    parser_events_stats = subparsers.add_parser('events_stats',
                                        help='NEW FEATURE -- get read stats of events.')

    parser_events_stats.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')

    parser_events_stats.add_argument('-f5', '--fast5type',
                              dest='f5type',
                              metavar='STRING',
                              choices=['minknow', 'metrichor'],
                              default='minknow',
                              help='''Specify the type of fast5 file - minknow output, or metrichor base-called output.
                                        choices = minknow, metrichor
                                        The events have different hdf5 locations based on type.
                                        Default: minknow
                                        Minknow columns: start, length, mean, variance
                                        Metrichor columns: mean, std dev, start, length''')
    parser_events_stats.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" where extension can be only pdf and jpg for now.''',
                             default=None)

    
    parser_events_stats.set_defaults(func=run_subtool)



    ##########
    # get_events
    ##########
    parser_get_events = subparsers.add_parser('get_events',
                                        help='NEW FEATURE -- get read stats of events.')
    parser_get_events.add_argument('-i', '--input',
                                   required=True,
                                   type=str,
                                   help='''Provide path to fast5 file.''')
    
    parser_get_events.add_argument('-f5', '--fast5type',
                              dest='f5type',
                              metavar='STRING',
                              choices=['minknow', 'metrichor'],
                              default='minknow',
                              help='''Specify the type of fast5 file - minknow output, or metrichor base-called output.
                                        choices = minknow, metrichor
                                        The events have different hdf5 locations based on type.
                                        Default: minknow
                                        Minknow columns: start, length, mean, variance
                                        Metrichor columns: mean, std dev, start, length''')
                                   
    
    parser_get_events.add_argument('-out', '--outputformat',
                             dest='out_fmt',
                             type=str,
                             help='''Return events in column structure of minknow or metrichor (defined in fast5type description).
                                    Default is output format of f5 type chosen.''',
                             default=None)
    parser_get_events.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" ''',
                             default=None)

    
    parser_get_events.set_defaults(func=run_subtool)


     ##########
    # get_model
    ##########
    parser_get_model = subparsers.add_parser('get_model',
                                        help='NEW FEATURE -- get read stats of events.')

    parser_get_model.add_argument('-i', '--input',
                                   required=True,
                                   type=str,
                                   help='''Provide path to fast5 file. Has to be metrichor output f5 file.''')

    parser_get_model.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" ''',
                             default=None)

    
    parser_get_model.set_defaults(func=run_subtool)

                                   
    

     ##########
    # get_metadata
    ##########
    parser_get_metadata = subparsers.add_parser('get_metadata',
                                        help='NEW FEATURE -- get read stats of events.')

    parser_get_metadata.add_argument('-i', '--input',
                                   required=True,
                                   type=str,
                                   help='''Provide path to fast5 file.''')
    
    parser_get_metadata.add_argument('-f5', '--fast5type', 
                              dest='f5type',
                              metavar='STRING',
                              choices=['minknow', 'metrichor'],
                              default='minknow',
                              help='''Specify the type of fast5 file - minknow output, or metrichor base-called output.
                                        choices = minknow, metrichor
                                        The events have different hdf5 locations based on type.
                                        Default: minknow
                                        Minknow columns: start, length, mean, variance
                                        Metrichor columns: mean, std dev, start, length''')  
##    parser_get_metadata.add_argument('-out', '--outputformat',
##                             dest='out_fmt',
##                             type=str,
##                             help='''Return events in column structure of minknow or metrichor (defined in fast5type description).
##                                    Default is output format of f5 type chosen.''',
##                             default=None)
    parser_get_metadata.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" ''',
                             default=None)

    
    parser_get_metadata.set_defaults(func=run_subtool)





    ##########
    # winner -- make poretools winner better for me.
    ##########
    parser_winner = subparsers.add_parser('winner',
                                        help='Get the longest read from a set of FAST5 files')
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
    # fragstats
    ##########
    parser_fragstats = subparsers.add_parser('fragstats',
                                        help='''Returns tab-delimited table with columns:
                                    1=readname, 2=estimated fragment size, 3=number input events, 4=if complement detected, 5=if 2D detected,
                                    6=num template events, 7=num complement events, 8=length of 2D sequence, 9=length of template sequence,
                                    10=length of complement sequence, 11=mean qscore of 2D sequence, 12=mean qscore of template sequence,
                                    13=mean qscore of complement.
                                    Estimates fragment size in the following way.
                                    If has 2D, this is length of 2D read.
                                    If template only, it is the length of template.
                                    If template and complement, but no 2D, it is length of the longer read between template and complement.
                                                    ''')
    parser_fragstats.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_fragstats.set_defaults(func=run_subtool)

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
    # get num events
    ##########
    parser_numevents = subparsers.add_parser('numevents',
                                        help='Print 2 column list of file and number of input events in file.')
    parser_numevents.add_argument('files', metavar='FILES', nargs='+',
                               help='The input FAST5 files.')
    parser_numevents.set_defaults(func=run_subtool)


    ##########
    # alignment
    ##########
    parser_align = subparsers.add_parser('align',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Performs alignments -- returns alignments, stats, plots')
    align_subparsers = parser_align.add_subparsers(title='[align-commands]', dest='align_command', parser_class=ArgumentParserWithDefaults)
    parser_align_fitting = align_subparsers.add_parser('fitting', help="fitting of 1 DNA seq to another")
    parser_align_blasr = align_subparsers.add_parser('blasr', help="BLASR")
    parser_align_blastn = align_subparsers.add_parser('blastn', help="BLASTN")
    parser_align_last = align_subparsers.add_parser('last', help="LAST")
    ## BLASR
    parser_align_blasr.add_argument("--plot")
    parser_align_blasr.add_argument("--blasr_file", type=str,default=None)
##    parser_align_filetype = parser_align.add_mutually_exclusive_group(required=True)
##    parser_align_filetype.add_argument('-fa', '--fasta',
##                              dest='fasta',
##                              default=None,
##                              type=str,
##                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
##    parser_align_filetype.add_argument('-fq', '--fastq',
##                              dest='fastq',
##                              default=None,
##                              type=str,
##                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_align.add_argument('--sequence',
                             dest='sequence',
                             type=str,
                             help='''Provide a sequence < min read length. Default is the KanRR fragment: CGTACGCTGCAGGTCG''',
                             default='CGTACGCTGCAGGTCG') ##shortest version of KanR fragment that maximized scores on some pilot reads
    parser_align.add_argument('-m', '--multiple-sequences',
                             dest='multiple_sequences',
                             type=str, default=None,
                             help='''Provide a path to a file with 1 sequence per line.
                                    For each read in the fastx file, it will report the fitting alignment for the
                                    sequence in this file with the best fitting aln.
                                    Each time it encounters a score that ties the current max score, it exchanges the older fiting aln
                                    info for the new fitting aln info with a 50%% probability.
                                    This way there is a random assignment of the best barcode.
                                    Use --all-scores instead to get an output with all max scores and barcodes returned.''')
    parser_align.add_argument('-w', '--with-read-names',
                             dest='with_read_names', action="store_true",
                             default=False,
                             help='''If set, will print "readname, startPosInRead, fitAlnScore, fitAlnScore/queryLen";
                                else just "startPosInRead,fitAlnScore, fitAlnScore/queryLen".
                                Start position is in pythonese (0-based).''')
    parser_align.add_argument('-e', '--with-edit-distances',
                             dest='with_edit_distances', action="store_true",
                             default=False,
                             help='''If set, edit dist will be incl in output''')
    parser_align.add_argument('-a', '--with-aln-seqs',
                             dest='with_aln_seqs', action="store_true",
                             default=False,
                             help='''If set, the aligned versions of sequences 1 (read) and 2 (provided) will be printed.''')
    parser_align.add_argument('-r', '--random-sequence',
                             dest='random_sequence', type=int,
                             default=False,
                             help='''Provide integer for random sequence length. This option overrides --sequence.''')
    parser_align_seqtransform = parser_align.add_mutually_exclusive_group()
    parser_align_seqtransform.add_argument("-c", "--complement", action="store_true", default=False,
                                                 help=''' Use complement of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> TTGG''')
    parser_align_seqtransform.add_argument("-rc", "--reverse_complement", action="store_true", default=False,
                                                 help=''' Use reverse complement of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> GGTT''')
    parser_align_seqtransform.add_argument("-rs", "--reverse_sequence", action="store_true", default=False,
                                                 help=''' Use reverse sequence of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> CCAA''')

    parser_align.set_defaults(func=run_subtool)


                                     


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
