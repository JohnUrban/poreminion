poreminion - 0.4.2
==========

Additional tools for analyzing Oxford Nanopore minION data by John Urban. 

Tools:

    uncalled            Find Fast5 files that were not base-called.
    timetest            Find Fast5 files that have event times that are earlier than event times before it suggesting malfunction/erroneous read.
    fragstats           Run this on set of base-called fast5 files.
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
                        19 = num skips in template (is actually number 0 moves found in extensive analysis)
                        20 = num skips in complement (is actually number 0 moves found in extensive analysis)
                        21 = num stays in template (is actually number 2 moves found in extensive analysis, any 3,4,5 moves not counted here)
                        22 = num stays in complement (is actually number 2 moves found in extensive analysis, any 3,4,5 moves not counted here)
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
                        
                        If -g4/--quadruplex used:
                        Final+1 = number of G4 motifs in 2D read: '([gG]{3,}\w{1,7}){3,}[gG]{3,}' 
                        Final+2 = number of G4 motifs in template read 
                        Final+3 = number of G4 motifs in complement read
                        Final+4 = number of G4 complement motifs in 2D reads: '([cC]{3,}\w{1,7}){3,}[cC]{3,}'
                        Final+5 = number of G4 complement motifs in template read (i.e. inferred complement strand count given template read)
                        Final+6 = number of G4 complement motifs in complement read (i.e. inferred template strand count given complement read)
                        
                        If --checktime used:
                        Final column (after even G4 info) = 0 or 1 for no/yes there is a time error present.
                        
                        
                        Estimates molecule/fragment size in the following way.
                        If has 2D, molecule size is the length of 2D read.
                        If template only, molecule size is the length of template read.
                        If template and complement, but no 2D, molecule size is length of the longer read between template and complement.
                        Molecule size allows calculation of total non-redundant data.
                        This is the sum of unique molecule lengths rather than summing all read types from each molecule.
                        From the molecule sizes, the "Molecule N50" can be computed using the nx subcommand on the fragstats file and specifying colum 2.
                                                                            
    fragsummary         To summarize fragstats, use this with a tab-delimited, fragstats table file (output of fragstats subcommand).
    nx                  Computes N50 or NX values on columns of a file or from comma-separated list.
    fragrobust          Looks at fragsizes in fragstats. Sees what percent of fragsizes are "robust" to all sequence lengths from same molecule.
    pct2d               Get the proportion of reads that have a 2D read
    has2d               Prints 2 columns: filename, has2D =  True/False
    numevents           Print 2 column list of file and number of input events in file.
    events              Look at events inside raw and basecalled fast5 files. 
    staypos             Get BED output of stay positions in read(s). 
    info                Get info about run and, if file is basecalled, basecalling. 
    g4                  Use quadparser suite (for identifying G4 motifs) on set of fast5 files (or in a FASTA/FASTQ file) and get a BED file with info for each match.
                            The default parameters search for '([gG]{3,}\w{1,7}){3,}[gG]{3,}' and its complement '([cC]{3,}\w{1,7}){3,}[cC]{3,}'.
                            See: http://en.wikipedia.org/wiki/G-quadruplex#Quadruplex_prediction_techniques
                            
                            This automates the regex sub-command to search for G4s with given paramters.
                            See regex for more info on output and searching for any regular expression.
                                
    regex               Search sequences in set of fast5 files (or in a FASTA/FASTQ file) for a regular expression.
                            Output BED file has default columns:
                            1. Name of sequence 
                        
                            2. Start of the match 
                        
                            3. End of the match
                            4. Strand (+/- relative to sequence given, NOT to be confised with template/complement reads.)
                            5. Optional Matched sequence (--reportseq/-s)
                        
                            These can be changed with --outformat/-o which allows you to report name,start,end,strand,seq in any order.
                        
                            If --counts is used, default columns are:
                            1. name
                            2. pos strand count
                            3. neg strand count
                            4. total count
                            
                            This script will write out all positive strand entries of a given sequence followed by all negative strand entries.
                            If name,start,end are used as first 3 columns, sortBed from BEDtools (or unix sort) can sort the BED file based on coordinates if needed.
                            
    dataconc            Plot sum of read lengths in each bin for a given set of bins for a set of FAST5 files.
                        This is the type of plot seen in MinKNOW while sequencing.
    qualpos             Get the qual score distribution over positions in reads
    qualdist            Get the qual score composition of a set of FAST5 files.
                        This tool is from poretools, but poreminion allows you to select the type of read.
    kmer                Count kmers in reads or reference.
    kmerplot            Plot kmer counts in reads or reference.
    kmerdiff            Get fold-enrichment values of kmers in reads vs reference.
    winner              Get the longest read from a set of FAST5 files.
                        Similar to poretools winner, only allows type=each and offers a details only option.
    seqlen              Get sequence lengths from set of FAST5 files.
                        By default it will attempt to give read lengths for template, complement, and 2d.
                        Use optional flags to exclude any of these read types.


Requirements
==========

poretools (https://github.com/arq5x/poretools)

HDF5 >= 1.8.7 (http://www.hdfgroup.org/HDF5/)

R >= 3.0.0

Python >= 2.7

rpy2 >= 2.4.2

h5py >= 2.0

pandas >= 0.14.1

matplotlib >= 1.4.0

edgeR (only for kmer 'differential expression' analysis)

Note
======
PoreMinion contains some tools that have been made on top of Aaron Quinlan's and Nick Loman's poretools and others that use other software.
PoreMinion will requiring additional dependencies on top of what poretools requires (specifically pandas and matplotlib libraries). 


INSTALL:
=======
download zip

cd porminion-master/

python setup.py install


Dependency installation (assuming Mac OS X):
-------------------------------------------
brew install hdf5

pip install h5py

download poretools zip; cd poretools; python setup.py install

pip install pandas

pip install matplotlib

