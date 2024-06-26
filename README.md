# Taiji2
Major update of Taiji with enhanced downstream analysis, featuring:

- TF’s Edgeweight per locus (I think this is really powerful to understand the neighboring TF effect and visualization)


- TF wave


- TF-regulatee analysis
  - Script for TF-regulatee edgeweight
  - comparative TF-regulatee landscape and clustering (Ext. Fig. 9 related)

- TF-TF interaction visualization. 
  - Script for Fig. 6b.


- TF community construction and visualization
Inspired by DBPNet72, which is a framework to identify cooperations between DNA-binding proteins using Chromatin immunoprecipitation followed by sequencing (ChIP-seq) and Hi-C data, we constructed the TF interaction network based on Taiji’s output, which is TF-regulatee network. For each context, we first combined the cell state-important TFs and cell state-specific TFs. In total, 159 TFs for TEXterm and 170 TFs for TRM were selected. 

- Perturb seq integration/heuristic score calculation (additional input required).
  -Input: TF KO DEGs + Taiji edgeweight
  -Output: Heuristic score


