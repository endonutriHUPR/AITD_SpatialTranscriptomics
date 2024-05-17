####### FROM SEURAT TO ANNDATA #######
  
  
  ##### PACKAGES
  library(Seurat)
  library(SeuratObject)
  Sys.setenv(RETICULATE_PYTHON="/home/endoprincesa/anaconda3/envs/sceasy/bin/python3")
  library(reticulate)
  reticulate::use_condaenv("/home/endoprincesa/anaconda3/envs/sceasy/")
  library(sceasy)
  library(dplyr)

##### DATA
together_sc <- readRDS("sc_analysis.rds")

  
setdiff (rownames(together_sc@meta.data), colnames (together_sc@assays$RNA@counts))

together_sc@meta.data <- together_sc@meta.data %>% filter(!(row.names(together_sc@meta.data) %in% c(setdiff(rownames(together_sc@meta.data), colnames(together_sc@assays$RNA@counts)))) ) 

sceasy::convertFormat(together_sc, from="seurat", to="anndata",main_layer = "counts",
                      outFile='/sc_analysis_new.h5ad')
