poreminion
==========

Additional tools for analyzing Oxford Nanopore minION data (built on top of poretools)

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
poreminion contains tools that have been made with Aaron Quinlan and Nick Loman's poretools in mind, but with requiring additional dependencies on top of what poretools requires (specifically pandas and matplotlib libraries). Some or all of the additional tools in poreminion are also available in my fork of poretools and may eventually become integrated into the main poretools if the authors are interested. For now, poreminion will contain strictly additional functionality on top of poretools (without overlapping features/analyses for nanopore reads).

What's in a name?
=================
poreminion is named such because it is a minion to (dependent on) poretools and because it works with minION data. Therefore it is a play on Oxford Nanopore's play on the word minion. Super meta!
