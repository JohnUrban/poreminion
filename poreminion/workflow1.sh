#!/bin/bash

## provide F5 Dir as $1, and genome/ref seq for aln analysis for $2
## provide a reference as $2 -- for Sciara this might be PacBio contigs for now


##stats
echo STATS
for type in all 2D fwd rev; do
poretools stats --type $type $1 > stats.${type}.txt
done

## hist
echo HIST
poretools hist $1 --saveas hist.png

## yield
echo YIELD
poretools yield_plot $1 --saveas yield.png

## data conc as percent
echo DATA_CONC
poreminion data_conc --percent $1 --saveas dataconc_percent.jpg
## data conc as cumulative percent
poreminion data_conc --percent --cumulative $1 --saveas dataconc_percent.jpg


## qual pos by type
echo QUAL_POS
for type in all 2D fwd rev; do
poreminion qualpos --type $type $1 --saveas readquality.${type}.jpg
done

## winner sizes -- depends on faSize for now -- can just make poreminion do this
## winner fasta
echo WINNERS
poreminion winner --type each --details ../../fast5/downloads/filtered/ > winners.txt
poreminion winner --type each --details ../../fast5/downloads/filtered/ > winners.fa

## for type in all 2D fwd rev; do poretools winner --type $type $1 | faSize -detailed stdin >> winners.txt; done
## for type in 2D fwd rev; do poretools winner --type $type $1 > winner.${type}.fa; done

## winner alns with blasr
echo WINNER ALN TO REF
blasr winners.fa $2 -out winners.blasr -bestn 1
awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8-$7,100*($8-$7)/$9,$10,$11,$12,$11-$10,100*($11-$10)/$12}' winners.blasr >> winners.blasr.tmp
mv winners.blasr.tmp winners.blasr



## OLD CODE TO DELETE
#for type in 2D fwd rev; do
#  blasr winner.${type}.fa $2 -out winner.${type}.blasr
#  echo "#qname tname qStrand tStrand Score percent_similarity tStart tEnd tLength tAlnLength %t_cov qStart qEnd qLength qAlnLength %q_cov" > winner.${type}.blasr.tmp
#  awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8-$7,100*($8-$7)/$9,$10,$11,$12,$11-$10,100*($11-$10)/$12}' winner.${type}.blasr >> winner.${type}.blasr.tmp
#  mv winner.${type}.blasr.tmp winner.${type}.blasr
#done


## Other alignment stuff....
#COMING SOON


## KMER ANALYSIS Vs. Reference
echo KMER ANALYSIS
poreminion kmer $1 --saveas 5mer.reads
poreminion kmer $1 --fasta $2 --rev-comp --saveas 5mer.rc.reference
## plot ref (x-axis) vs reads (y-axis)
poreminion kmerplot -t1 5mer.rc.reference -t2 5mer.reads --saveas 5mer.scatter.jpg
## diff abundance, volcano, smear, ntcontent of diff abund kmers
poreminion kmerdiff -t1 5mer.rc.reference -t2 5mer.reads --saveas 5mer.diffabundance.txt --volcano 5mer.volcano.jpg --smear 5mer.smear.jpg --nt-content fc:2 > 5mer.ntcontent.txt
## top 10 enrched 5mers
sort -k 3,3nr 5mer.diffabundance.txt | head | cut -f 2,3,4,5,6 > 5mer.top10enrichedInReads.txt
sort -k 3,3n 5mer.diffabundance.txt | head | cut -f 2,3,4,5,6 > 5mer.top10depletedInReads.txt







##BLASR -- https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
## default output -- qName tName qStrand tStrand score percentSimilarity tStart tEnd tLength qStart qEnd qLength nCells
# space-delim -- 11 fields
# blasr winner.fa genome.fa  
# blasr winner.fa genome.fa -m 1

##SAM output
# blasr winner.fa genome.fa  -sam 
#  (1) "XS": read alignment start position without counting previous soft clips 
#  (2) "XE": read alignment end position without counting previous soft clips 
#  (3) "XL": aligned read length
#  (4) "XQ": query sequence length
#  (5) "XT": number of continues reads, always 1 for blasr


## human readable aln fmt
# blasr winner.fa genome.fa  -m 0 

## XML fmt
# blasr winner.fa genome.fa  -m 2

## vulgar fmt --- represents sequence as symbols for match, mismatch, gap, insertion, del, etc -- also has info in beginning
# blasr winner.fa genome.fa  -m 3

## 13 field space delim -- qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
# blasr winner.fa genome.fa  -m 4 

## 19 field space delim -- qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
# blasr winner.fa genome.fa  -m 5
