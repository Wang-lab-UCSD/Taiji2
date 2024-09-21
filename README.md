# Taiji2

Major update of Taiji with enhanced downstream analysis, featuring: 

-   Cell-state specificity analysis

-   TF’s edgeweight per locus calculation

-   TF transcriptional wave

-   TF-regulatee analysis

-   TF-TF interaction network

-   Perturb seq integration/heuristic score calculation.


<img src="https://github.com/cong-003/Taiji2/blob/main/figures/summary_fig.png" width="800"/>

## Instructions

### run Taiji pipeline for paired RNA-seq and ATAC-seq data

First install Taiji. Check [Taiji github](https://taiji-pipeline.github.io/)

``` bash
curl -L https://github.com/Taiji-pipeline/Taiji/releases/latest/download/taiji-CentOS-x86_64 -o taiji
chmod +x taiji
./taiji --help
```

Then prepare configure file and input file following instructions in [Taiji website](https://taiji-pipeline.github.io/)

``` bash
taiji run --config config.yml -n 3 +RTS -N3
```

To replicate the paper's results, use the configure file and input file in this [repo](https://github.com/Wang-lab-UCSD/Taiji2/tree/main/inputs/).

### run downstream analysis

After running Taiji, you will have `GeneRanks.tsv` file in `/some_path_to_Taiji/output/`, which stores the PageRank scores of TFs across samples. This will be the major data in the following downstream analysis.

Additionally, gene expression (normalized by TPM) file `expression_profile.tsv` is also available in folder `/some_path_to_Taiji/output/RNASeq/`

#### cell state specificity analysis

In addition to PageRank scores and gene expression matrices, group file with cell state annotation is also needed. Check the [demo group file](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig1/group_file.txt) for example.

Follow the [tutorial](https://rpubs.com/cong003/1201450) step by step. Required R packages are listed in the beginning. Alternatively, you can download the original [code](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig1/Fig1_TF_activity_visual.Rmd) for further customization.

#### TF transcriptional wave

TF transcriptional waves are different patterns patterns of TF groups which govern different differentiation pathways. To construct waves, you need to have pre-defined layout for differentiation path. The example dataframe is shown below, where `x` and `y` are coordinates of each cell state; `labelposx` and `labelposy` are label coordinates.

``` r
wavedf <- data.frame(x = c(1,2,2,3,3,4,4,5,5), 
                     y = c(3,5,1,3,1,3,0,4,2),
                     samplename = c("Naive","TE","TexProg","MP","TexInt","TRM","TexTerm","TEM","TCM"),
                     labelposx = c(1,2,2,3,3,4,4,5,5),
                     labelposy = c(3,5,1,3,1.4,3,0,4,2)-0.4)
```

Follow the [tutorial](https://rpubs.com/cong003/1201435) step by step. Alternatively, you can download the original [code](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig1/Fig1_TF_wave.Rmd) for further customization.

#### TF-regulatee network analysis

Taiji generated the regulatory network showing the regulatory relationship between TF and its target genes (regulatee) with edge weight, which represents the regulatory strength. The TF-regulatee network files "/some_path_to_taiji/output/Network/sample_name/edges_combined.csv". Each sample has its own network file.

Based on this, we can calculate the log2 fold change of edge weights between two cell states, in our case, between TexTerm and TRM. See the [example file](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/log2FC_btw_TEX_and_TRM_mean_edge_weight_subset_TFs_v2_top500_genes.csv)

Follow the [tutorial](https://rpubs.com/cong003/1201748) step by step. Alternatively, you can download the original [code](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/differential_edge_weight_heatmap.Rmd) for further customization.

#### TF-TF interaction network

First, we constructed the TF communities. From the above TF-regulatee network, we can derive TF-TF correlation based on the edge weight profile. Check the file format [here](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/TF_corr_TRM_subset_TFs_v2.csv).

Follow the [tutorial](https://rpubs.com/cong003/1216208) step by step. Alternatively, you can download the original [code](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/community_construction.Rmd) for further customization.

Next, we visualized the TF communities following the [tutorial](https://rpubs.com/cong003/1201736). The original code can be downloaded [here](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/TF_TF_network_visual.Rmd)

Functional analysis can be performed as well. Check the [tutorial](https://rpubs.com/cong003/1216075) and raw [code](https://github.com/Wang-lab-UCSD/Taiji2/blob/main/figures/Fig6/community_analysis.Rmd)


#### perturb-seq analysis
    
## Resources

-   [paer website]()
-   [Taiji website](https://taiji-pipeline.github.io/)
-   download Seurat objects [here](https://ucsdcloud-my.sharepoint.com/:f:/g/personal/ajambor_ucsd_edu/Eh-PQxt5WxJHjo5whw01KqYB1vOc-BFTlutg2Var8xzfeQ?e=TbBrpE) and move to appropriate folders to run tutorials

## How to cite

Chung, H. Kay, Cong Liu, Ming Sun, Eduardo Casillas, Timothy Chen, Brent Chick, Jun Wang et al. "Multiomics atlas-assisted discovery of transcription factors enables specific cell state programming." bioRxiv (2023): 2023-01. <https://www.biorxiv.org/content/10.1101/2023.01.03.522354v3.abstract>.
