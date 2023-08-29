# Autoimmune Thyroid Diseases with Spatial Transcriptomics

This repository host the code used in our paper titled "Epithelial and stromal landscape heterogeneity at spatial resolution in autoimmune thyroid diseases". 

## Abstract

Autoimmune thyroid diseases (AITD) such as Graves' disease (GD) or Hashimoto's thyroiditis (HT) are autoimmune, organ-specific diseases related to a complex and multifactorial interplay of specific susceptibility genes and environmental exposures. T lymphocytes and their secretory cytokines  play indispensable roles in disrupting tolerance. However, their roles are often complex and full of interactions among distinct components of the thyroid ecosystem. Novel technologies such as spatial transcriptomics (ST), allow us to  explore the molecular architecture and the heterogeneity and landscape of the different cell-states within tissues.  Here, we apply ST to focus on changes within the stroma and epithelial cells involved in these pathologies to further understand the complexity of AITD. We profile 8 human thyroid (3 HT, 3 GD and 2 controls) by ST and spatially localized and annotated  TFCs, connective tissue and vessel areas to be independently studied. We identified damaged antigen-presenting TFCs located close to the immune infiltration that contributed to the autoperpetuation of the immune response in AITD. Interestingly, the combination of single cell data and ST revealed the participation of different fibroblast subpopulations in HT and GD patients. Specifically, GD-associated ADIRF+ myofibroblasts surrounding TFCs and inflammation-associated fibroblasts (IAFs) within the connective tissue of HT patients. Furthermore, the presence of fenestrated PLVAP+ vessels in AITD, especially in GD, was highlighted. Our data provides the molecular and cellular heterogeneity of AITD microenvironment and points toward the importance of different novel markers associated to cell subpopulations in the stroma and thyroid epithelium that could be essential to understand the pathogenesis of AITD.

## Quality control and analysis

We can mainly summarize the bioinformatic protocol in four steps:

### 1) Sequencing quality control

_QC_seq_ST.R_ we overlooked the quality of the sequencing per sample. We got as output some tables and graphs for a easier evaluation.

### 2) Spatial quality control

Spots-based sequencing allows as to evaluate their quality in a spatial context. As an example, we can observe lower % of ribosome in TFCs areas compared to immune infiltration. _QC_spatial_ST.R_ was used for that aim.

### 3) Samples integration, clustering, differential expression (DE) analysis and enrichment analysis

We used the code available in the R script _integration_of_samples_and_analysis.R_.

### 4) Re-clustering

As developed in _reclustering.R_ , we firstly isolated the interesting cluster to 1) separate samples, normalize, dimmensional reduced them and 2) merge, scale and re-integrate them. Finally; UMAP, re-clustering and DE and enrichment analysis were performed.
