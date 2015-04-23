poreminion
==========

Additional tools for analyzing Oxford Nanopore minION data (built on top of poretools) by John Urban.

Tools:

    uncalled            Find Fast5 files that were not base-called.
    timetest            Find Fast5 files that have event times that are earlier than event times before it suggesting malfunction/erroneous read.
    fragstats           Run this on set of base-called fast5 files.
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
                                                                            
    nx                  Computes N50 or NX values on columns of a file or from comma-separated list.
    pct2d               Get the proportion of reads that have a 2D read
    has2d               Prints 2 columns: filename, has2D =  True/False
    numevents           Print 2 column list of file and number of input events in file.
    events              Look at events inside raw and basecalled fast5 files. 
    dataconc            Plot sum of read lengths in each bin for a given set of bins for a set of FAST5 files.
                        This is the type of plot seen in MinKNOW while sequencing.
    qualpos             Get the qual score distribution over positions in reads
    qualdist            Get the qual score composition of a set of FAST5 files
    kmer                Count kmers in reads or reference.
    kmerplot            Plot kmer counts in reads or reference.
    kmerdiff            Get fold-enrichment values of kmers in reads vs reference.
    winner              Get the longest read from a set of FAST5 files.
                        Similar to poretools winner, only allows type=each and offers a details only option.

Requirements
==========

poretools (https://github.com/arq5x/poretools)

HDF5 >= 1.8.7 (http://www.hdfgroup.org/HDF5/)

R >= 3.0.0

Python >= 2.7

rpy2 >= 2.4.2

h5py >= 2.0

pandas>=0.14.1

matplotlib>=1.4.0

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
