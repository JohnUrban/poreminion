poreminion
==========

Additional tools for analyzing Oxford Nanopore minION data (built on top of poretools) by John Urban.
Warning: I only work on PoreMinion every once in a while and some tools are un-finished and/or broken.
Current tools I know to work are:
	- data_conc
	- qualpos
	- qualdist
	- kmer
	- kmerplot 
	- kmerdiff
	- winner
	- pct2d
	- uncalled
	- numevents
The 'align' subtool is something I tooled with, which was just going to be a wrapper over various aligners. I believe it is not in working order.
Others such as events_stats, get_events, get_model, get_metadata are also un-finished. I have similar functions in another suite that I will make available soon.

- data concentration (DC) plots:      
-         poreminion data_conc

- cumulative data concentration plots:    
-         poreminion data_conc --cumulative

- DC plots as % of data:
-         poreminion data_conc --percent  
-         poreminion data_conc --cumulative --percent

- quality vs. position boxplots:
-         poreminion qualpos

-   of only 2D reads: 
-         poreminion qualpos --type 2D

-   of only template reads: 
-         poreminion qualpos --type fwd

-   of only template reads: 
-         poreminion qualpos --type rev


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
