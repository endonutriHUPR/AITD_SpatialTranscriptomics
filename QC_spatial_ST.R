# ============================================================================= #
# ============================================================================= #
########################### QC spatial metric summary ###########################
# ============================================================================= #
# ============================================================================= #

# ============================================================================= #
############################### 0. Initialization ###############################
# ============================================================================= #

########### A) Set working directory

setwd("~/Spatial_gene_expression/")

# New working directory
dir.create("./QC_spatial_level_G2_C2")
setwd("./QC_spatial_level_G2_C2")


########### B) Libraries

library(Seurat)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)

########### C) Sample IDS

# EXAMPLE #
###!!!### Make sure about the order

samples_ids <- c("GD2", "C2", "GD3", "C1")

########### D) Directories where the files are (4)

GD2 <- Seurat::Load10X_Spatial(data.dir = "~/AITD/Spatial_gene_expression/Mapping/G2_C2/GD2/outs/",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = samples_ids[1],
                                  filter.matrix = TRUE)

C2 <- Seurat::Load10X_Spatial(data.dir = "~/AITD/Spatial_gene_expression/Mapping/G2_C2/C2/outs/",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = samples_ids[2],
                                   filter.matrix = TRUE)

GD3 <- Seurat::Load10X_Spatial(data.dir = "~/AITD/Spatial_gene_expression/Mapping/G2_C2/GD3/outs/",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = samples_ids[3],
                                   filter.matrix = TRUE)

C1 <- Seurat::Load10X_Spatial(data.dir = "~/AITD/Spatial_gene_expression/Mapping/G2_C2/C1/outs/",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = samples_ids[4],
                                   filter.matrix = TRUE)



# ============================================================================= #
############################ 1. Some transformations ############################
# ============================================================================= #

########## Adding sample name

GD2[["sample_id"]] <- samples_ids[1]
C2[["sample_id"]] <- samples_ids[2]
GD3[["sample_id"]] <- samples_ids[3]
C1[["sample_id"]] <- samples_ids[4]


########## Combining the four objects
  
  se_obj <- merge(GD2, y = list(C2,GD3,C1),
                  add.cell.ids = samples_ids,
                  project = "BATLLE")
  
  ########## Add mitochondrial and ribosomal %
  
  se_obj[["percent.mito"]] <- Seurat::PercentageFeatureSet(
      object = se_obj,
      pattern = "^MT-")
    
  se_obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
      object = se_obj,
      pattern = "^RPL|^RPS")
  
  
  
  ########## Remove empty genes¨
  
  
  keep_genes <- rowSums(as.matrix(se_obj@assays$Spatial@counts)) != 0
  se_obj <- se_obj[keep_genes, ]
  
  
  
  
  # ============================================================================= #
  ############################## 2. Basic features  ###############################
  # ============================================================================= #
  
  ########### Unique genes per spot 
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]], 
                            ggplot2::aes(nFeature_Spatial),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::ggtitle("Unique genes per spot") +
    ggplot2::labs(x = "Number of Detected Genes",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  png("./unique_genes_per_spot.png", width = 700, height = 700)
  p1
  dev.off()
  
  ########## Plot number of features
  
  p1 <- lapply(samples_ids, function(id) {
      Seurat::SpatialFeaturePlot(
      object = se_obj[, se_obj$sample_id == id],
      features = "nFeature_Spatial",
      images = id) +
      ggplot2::labs(title = id) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18))
    
  }) 
    
  png("./nFeatures_plot.png", width = 700, height = 700)
  cowplot::plot_grid(
      plotlist = p1,
      align = "hv",
      axis = "trbl")
  
  dev.off()
  
  ######### Number of reads
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]],
                            ggplot2::aes(nCount_Spatial),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::ggtitle("Total counts per spots") +
    ggplot2::labs(x = "Library Size (total UMI)",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  png("./counts_per_spot.png", width = 700, height = 700)
  p1
  dev.off()
  
  ######### Plot spatial counts
  p1 <- lapply(samples_ids, function(id) {
    Seurat::SpatialFeaturePlot(
      object = se_obj[, se_obj$sample_id == id],
      features = "nCount_Spatial",
      images = id) +
      ggplot2::labs(title = id) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18))
    
  }) 
  
  png("./nCounts_plot.png", width = 700, height = 700)
  cowplot::plot_grid(
    plotlist = p1,
    align = "hv",
    axis = "trbl")
  
  dev.off()
  
  
  ######### Counts per gene
  
  count_mtrx <- Seurat::GetAssayData(object = se_obj,
                                     slot = "counts",
                                     assay = "Spatial")
  
  gene_attr <- lapply(samples_ids, function(id) {
    mask <- stringr::str_detect(
      string = colnames(count_mtrx),
      pattern = id)
    
    gene_attr <- data.frame(
      nUMI = Matrix::rowSums(count_mtrx[, mask]), 
      nSpots = Matrix::rowSums(count_mtrx[, mask] > 0),
      #gem_id = gem_id,
      sample_id = id)}) %>% dplyr::bind_rows()
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = gene_attr,
                            ggplot2::aes(nUMI),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle("Total counts per gene (log10 scale)") +
    ggpubr::theme_pubr()
  
  png("./counts_per_gene.png", width = 700, height = 700)
  p1
  dev.off()
  
  ####### Genes distributed in spots
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = gene_attr,
                            ggplot2::aes(nSpots),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::ggtitle("Total spots per gene") +
    ggpubr::theme_pubr()
  
  png("./genes_in spots.png", width = 700, height = 700)
  p1
  dev.off()
  
  ######### % Mitochondrial
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]], 
                            ggplot2::aes(percent.mito),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::ggtitle("Mitochondrial % per spot") +
    ggplot2::labs(x = "Mitochondrial % ",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  png("./percent_mito.png",width = 700, height = 700 )
  p1
  dev.off()
  
  # PLOT
  
  p1 <- lapply(samples_ids, function(id) {
    Seurat::SpatialFeaturePlot(
      object = se_obj[, se_obj$sample_id == id],
      features = "percent.mito",
      images = id) +
      ggplot2::labs(title = id) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18))
    
  })
  
  png("./percent_mito_spatial.png", width = 700, height = 700)
  cowplot::plot_grid(plotlist = p1)
  dev.off()
  
  
  ######### % Ribosomal
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]], 
                            ggplot2::aes(percent.ribo),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::facet_wrap(. ~ sample_id, scales = "free") +
    ggplot2::ggtitle("Ribosomal % per spot") +
    ggplot2::labs(x = "Ribosomal % ",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  png("./percent_ribo.png",width = 700, height = 700 )
  p1
  dev.off()
  
  # PLOT
  
  p1 <- lapply(samples_ids, function(id) {
    Seurat::SpatialFeaturePlot(
      object = se_obj[, se_obj$sample_id == id],
      features = "percent.ribo",
      images = id) +
      ggplot2::labs(title = id) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18))
    
  })
  
  png("./percent_ribo_spatial.png", width = 700, height = 700)
  cowplot::plot_grid(plotlist = p1)
  dev.off()
  
  
  
  # ============================================================================= #
  ############################## 3. Save the object  ##############################
  # ============================================================================= #
  
  lapply(samples_ids, function(id) {
    print(id)
    se_obj_sub <- se_obj[, se_obj$sample_id == id]
    sample_id <- unique(se_obj_sub$sample_id)
    
    # Remove other images
    se_obj_sub@images <- se_obj_sub@images[Seurat::Images(se_obj_sub) == id]
    
    saveRDS(
      object = se_obj_sub,
      #file = here::here(glue::glue("./qc_se_{id}.rds")))
      file = paste0("./qc_se_", id, ".rds"))
  })
    
