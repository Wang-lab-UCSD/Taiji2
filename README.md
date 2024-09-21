# Taiji2
Major update of Taiji with enhanced downstream analysis, featuring:

- TF’s edgeweight per locus calculation

- TF-regulatee analysis
  - [differential edge weight analysis](https://rpubs.com/cong003/1201748)
 
- TF-TF interaction network
  - [TF community construction](https://rpubs.com/cong003/1216208)
  - [community visualization](https://rpubs.com/cong003/1201736)
  - [functional analysis](https://rpubs.com/cong003/1216075)
    
- TF transcriptional wave (additional input required).
  - Input: pre-defined differentiation path + Taiji pagerank
  - Output: transcriptional wave patterns
  - [tutorial](https://rpubs.com/cong003/1201435)

- Perturb seq integration/heuristic score calculation (additional input required).
  - Input: TF KO DEGs + Taiji edgeweight
  - Output: Heuristic score
 
- [Cell-state specificity analysis](https://rpubs.com/cong003/1201450)

<img src="https://github.com/cong-003/Taiji2/blob/main/figures/summary_fig.png" width="800">

## Instructions
### run Taiji pipeline for paired RNA-seq and ATAC-seq data
First install Taiji. Check [Taiji github](https://taiji-pipeline.github.io/)
```bash
curl -L https://github.com/Taiji-pipeline/Taiji/releases/latest/download/taiji-CentOS-x86_64 -o taiji
chmod +x taiji
./taiji --help
```


## Resources
- website:
- tutorial and demo: 
- intro of Taiji: https://taiji-pipeline.github.io/
- download Seurat objects [here](https://ucsdcloud-my.sharepoint.com/:f:/g/personal/ajambor_ucsd_edu/Eh-PQxt5WxJHjo5whw01KqYB1vOc-BFTlutg2Var8xzfeQ?e=TbBrpE) and move to appropriate folders to run tutorials

## How to cite
Chung, H. Kay, Cong Liu, Ming Sun, Eduardo Casillas, Timothy Chen, Brent Chick, Jun Wang et al. "Multiomics atlas-assisted discovery of transcription factors enables specific cell state programming." bioRxiv (2023): 2023-01. https://www.biorxiv.org/content/10.1101/2023.01.03.522354v3.abstract.
 
