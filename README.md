# Autoimmune Thyroid Diseases with Spatial Transcriptomics

This repository host the code used in our paper titled "Unraveling the molecular architecture of autoimmune thyroid diseases at spatial resolution". 

## Abstract

Autoimmune thyroid diseases (AITD) such as Graves' disease (GD) or Hashimoto's thyroiditis (HT) are organ-specific diseases that involve complex interactions between distinct components of the thyroid ecosystem. We applied spatial transcriptomics (ST) to explore the molecular architecture, heterogeneity and location of different cells present in the thyroid tissue, including thyroid follicular cells (TFCs), stromal cells such as fibroblasts, endothelial cells, and thyroid infiltrating lymphocytes. We identified damaged antigen-presenting TFCs that had an upregulated CD74/MIF axis in thyroids from AITD patients. Furthermore, we discerned two main fibroblast subpopulations in the connective tissue including ADIRF+ myofibroblasts, mainly upregulated in GD, and inflammatory-associated fibroblasts (IAFs), mainly overexpressed in HT patients. We also demonstrated an increase of fenestrated PLVAP+ vessels in AITD, especially in GD. Our data unveil novel stromal and thyroid epithelial cell subpopulations that could play an essential role in the pathogenesis of AITD.

![Graphical abstract SHORT](https://github.com/endonutriHUPR/AITD_SpatialTranscriptomics/assets/114569590/b0dd194a-2250-4524-9c07-f83096a7387b)

## Quality control and analysis

All this code was run in:<br>
R version 4.0.3 (2020-10-10)<br>
Platform: x86_64-conda-linux-gnu (64-bit)<br>
Running under: Ubuntu 20.04.6 LTS<br>
Seurat version 4.0.3.

We can mainly summarize the bioinformatic protocol in four steps:

### 1) Sequencing quality control

_QC_seq_ST.R_ we overlooked the quality of the sequencing per sample. We got as output some tables and graphs for a easier evaluation.

### 2) Spatial quality control

Spots-based sequencing allows as to evaluate their quality in a spatial context. As an example, we can observe lower % of ribosome in TFCs areas compared to immune infiltration. _QC_spatial_ST.R_ was used for that aim.

### 3) Samples integration, clustering, differential expression (DE) analysis and enrichment analysis

We used the code available in the R script _integration_of_samples_and_analysis.R_.

### 4) Re-clustering

As developed in _reclustering.R_ , we firstly isolated the interesting cluster to 1) separate samples, normalize, dimmensional reduced them and 2) merge, scale and re-integrate them. Finally; UMAP, re-clustering and DE and enrichment analysis were performed. 

## 5) Further analysis of Visium data
### Spots deconvolution

Regarding Python, version 3.9.18 was used. 

We performed the deconvolution of the spots in HT Visium samples using HRA001684 repository. To this aim, we chose cell2location package written in Python. Firstly, we transform R dataset to anndata using R package sceasy (_from_seurat_to_anndata_sceasy.R_). Then, we followed cell2location pipeline with our samples as in the script _cell2loc_st_ht.py_.

### Interactome
For ligand-receptor interaction we used the clusters we got in the analysis of our Visium data to run CellChat in its spatial data modality. All the steps are specified in _cellchat_spatial.R_ file.
