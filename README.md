# sciprog
This repo contains my final projects for the Scientific Programming course at BCG.

## Python Project
#### Package for disease gene prediction based on an protein-protein interaction network

The software `proj.py` within the `pyProj` directory provides a python function `genePred` which takes:
- a protein-protein interaction network from an individual species (e.g., human) from the STRING database, using only physical interactions,
- a set of known disease-related “seed” (or reference) genes,
- (optionally) a set of “candidate” disease genes (if not specified, all non-seed genes will be taken as candidates)
and ranks all “candidate” genes according to their vicinity to the seed genes on the network.  

It returns:
- a `results.txt` file containing the ranked list of the candidate genes, including the scores of their vicinty to the seed genes;
- a `resultsDetailed.txt` file with more detailed information on the results.
- (optionally) a `resultBarplot.pdf` file that graphically compares the scores of the individual candidate genes
- (optionally) a `resultNetwork.pdf` file that shows how distant or close seed genes and candidate genes are within the PPI network.  
  


## R Project
#### Determine the mutation type for a set of single nucleotide variants in a genome

The R package `mutType` within the `RProj` directory provides an R function `mutType` which takes:
- a set of mutations (single nucleotide variants, SNVs) in VCF format,
- the corresponding reference genome (e.g., human genome hg38),
- a parameter “context_length” which is a positive, odd integer  

It determines for each mutation (only SNVs; other mutations like indels are ignored) the corresponding mutation type as follows:
The mutation type is `UP[REF>ALT]DOWN` where
- `REF>ALT` is the single nucleotide variant from `REF` base to `ALT` base, e.g., “`C>T`”
- `UP` is one or more upstream bases from the reference genome (depending on the user parameter “context_length”)
- `DOWN is one or more downstream bases from the reference genome (same user parameter)

Optionally, the `mutTypeTable` function can be used to summarize the results and produce a barplot of the mutation types frequencies.

