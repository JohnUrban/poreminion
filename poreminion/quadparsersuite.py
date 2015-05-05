## Poreminion Quadparsersuite - Apr 24, 2015 JU
## My quadparsersuite code for fasta and fastq files originally appeared as an independent program
## at: https://github.com/JohnUrban/LexoNSseq2015/tree/master/quadparsersuite
## quadparsersuite was an expansion/modification of dario beraldi's original "quadparser"
## Dario's script Copy/pasted from:
## http://bioinformatics-misc.googlecode.com/svn-history/r16/trunk/quadparser.py
## on April 9, 2013 (11:42PM)
## Also see: http://code.google.com/p/bioinformatics-misc/source/browse/trunk/quadparser.py?spec=svn16&r=16

import re
import sys
import string
from poretools.Fast5File import *
from info import get_basename
from cStringIO import StringIO

def parseSequence(seq, seq_name, fwd_re, rev_re, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    for m in re.finditer(fwd_re, seq):
        formatDict = {'name':seq_name, 'start':m.start(), 'end':m.end(), 'strand':'+', 'seq':m.group(0)}
        if count_gtract:
            formatDict['gtracts'] = get_regex_count(formatDict['seq'], re.compile("[gG]{"+str(count_gtract)+",}"))
        print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])
    if noreverse is False:
        for m in re.finditer(rev_re, seq):
            formatDict = {'name':seq_name, 'start':m.start(), 'end':m.end(), 'strand':'-', 'seq':m.group(0)}
            if count_gtract:
                formatDict['gtracts'] = get_regex_count(formatDict['seq'], re.compile("[cC]{"+str(count_gtract)+",}"))
            print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])

def countSequence(seq, seq_name, fwd_re, rev_re, outformat, count_gtract=False, noreverse=False):
    ##count_gtract is not useful for this function for now, but needed to add it to play nice with other functions
    poscount = 0
    negcount = 0
    for m in re.finditer(fwd_re, seq):
        poscount += 1
    if noreverse is False:
        for m in re.finditer(rev_re, seq):
            negcount += 1
    formatDict = {'name':seq_name, 'pos':poscount, 'neg':negcount}
    print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])
    
def parseFasta(seq_fh, re_f, re_r, function, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    seq = []
    line = (seq_fh.readline()).strip()
    seq_name = re.sub('^>', '', line)
    line = (seq_fh.readline()).strip()
    while True:
        while line.startswith('>') is False:
            seq.append(line)
            line= (seq_fh.readline()).strip()
            if line == '':
                break
        seq = ''.join(seq)

        function(seq=seq, seq_name=seq_name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        
        seq_name = re.sub('^>', '', line)
        seq= []
        line= (seq_fh.readline()).strip()
        if line == '':
            break


def parseFastq(seq_fh, re_f, re_r, function, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    line = (seq_fh.readline()).strip()
    seq_name = re.sub('^@', '', line)
    seq = (seq_fh.readline()).strip()
    while True:
        function(seq=seq, seq_name=seq_name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        seq_fh.readline() ## skip + line
        seq_fh.readline() ## skip base call quality line
        line = (seq_fh.readline()).strip()   
        seq_name = re.sub('^@', '', line)
        seq = (seq_fh.readline()).strip()
        if line == '':
            break


def parseFast5(fast5dir, seqtype, re_f, re_r, function, outformat, count_gtract=False, noreverse=False):
    for fast5 in Fast5FileSet(fast5dir):
        fas = fast5.get_fastas(seqtype)
        for fa in fas:
            if fa is None:			
                continue
            basename,filename = fa.name.strip().split()
            filename = filename.split(".")[0].split("/")[-1]
            if filename.endswith("strand"):
                filename = filename[:-6]
            basename = basename.split("_")[-1]
            name = filename + basename
            function(seq=fa.seq, seq_name=name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        fast5.close()



    ## TODO -- output = name numG4_2d numG4temp numG4comp
    ## make "fragcombine" that combines frag style files - returns non-redundant ordered output
    ## perhaps its worthwhile for frag files to have header at top or assoc header file
    ## this will allow a program like fragcombine to have no prior (hard coded) knowledge (e.g. specified by options like -e -t) of order
##    f5=fast5.hdf5file
##    name = fast5.filename.split("/")[-1]

 
def get_regex_count(seq, regex):
    ## regex is a re.compile(string) object
    return len(re.findall(regex, seq))

def get_regex_count_in_fast5_fasta(fast5fasta, regex):
    if fast5fasta is not None:
        return get_regex_count(fast5fasta.seq, regex)
    else:
        return "-"
    
def get_regex_counts_in_fast5(fast5, regex, regex2=None):
    ## regex and regex2 are re.compile(string) objects
    ## often regex2 is the complement of regex (not necessary though)
    twoD = fast5.fastas.get("twodirections")
    template = fast5.fastas.get("template")
    complement = fast5.fastas.get("complement")
    counts = []
    counts.append(get_regex_count_in_fast5_fasta(twoD, regex))
    counts.append(get_regex_count_in_fast5_fasta(template, regex))
    counts.append(get_regex_count_in_fast5_fasta(complement, regex))
    if regex2 is not None:
        counts.append(get_regex_count_in_fast5_fasta(twoD, regex2))
        counts.append(get_regex_count_in_fast5_fasta(template, regex2))
        counts.append(get_regex_count_in_fast5_fasta(complement, regex2))
    return counts
        
    
def get_g4_regex(minG=3, maxN=7):
    return '([gG]{' + str(minG) + ',}\w{1,' + str(maxN) + '}){3,}[gG]{' + str(minG) + ',}'

def get_g4_revregex(minG=3, maxN=7):
    return '([cC]{' + str(minG) + ',}\w{1,' + str(maxN) + '}){3,}[cC]{' + str(minG) + ',}'

def get_complement_regex(regex, transtab=string.maketrans('actguACTGU', 'tgacaTGACA')):
    return regex.translate(transtab)



## should these return them compiled?

def run(parser, args):
    ## simplify fastx
    if args.fasta:
        args.fastx = args.fasta
    elif args.fastq:
        args.fastx = args.fastq
    else:
        args.fastx = None
        
    ## subcommand g4?
    if args.command == "g4":
        args.regex = get_g4_regex(args.minG, args.maxN)
        args.regexrev = get_g4_revregex(args.minG, args.maxN)

    ## define regexrev as complement of regex if none given
    if args.regexrev is None: ## then use complement of regex
        args.regexrev = get_complement_regex(args.regex)

    ## compile fwd and rev reg expressions
    re_f = re.compile(args.regex)
    re_r = re.compile(args.regexrev)

    ## process automated reporting/output format requests
    if args.reportseq:
        args.outformat = 'name,start,end,strand,seq'                                                      
    if args.counts:
        args.outformat = 'name,pos,neg'
    if args.numtracts:
        args.numtracts = int(args.minG)
        args.outformat += ",gtracts"

    ## allow fasta and fastq files be piped in as stdin
    if args.fastx == '-':
        fastx_fh = sys.stdin
        if len(args.fastx) > 1:
            sys.exit('\nPoreminion g4/regex: Only one input file at a time can be processed:\n--fasta: %s\n' %(args.input))
    elif not args.fast5:
        fastx_fh = open(args.fastx)


    ## determine function to use
    if args.counts:
        function = countSequence
    else:
        function = parseSequence


    ## execute  
    if args.fasta:
        parseFasta(fastx_fh, re_f, re_r, function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)
    elif args.fastq:
        parseFastq(fastx_fh, re_f, re_r, function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)
    elif args.fast5:
        parseFast5(args.fast5, args.type, re_f, re_r, function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)


