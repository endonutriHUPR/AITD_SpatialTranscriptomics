# ================================================================================= #
# ================================================================================= #
######################### INTEGRATION OF DATA WITH HARMONY ##########################
# ================================================================================= #
# ================================================================================= #



# ================================================================================= #
###################################### ARGUMENTS ####################################
# ================================================================================= #


#HT
HT1 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_H3_G1/qc_se_HT1.rds")
HT2 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_H3_G1/qc_se_HT2.rds")
HT3 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_H3_G1/qc_se_HT3_mod.rds")

#GD
GD1 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_H3_G1/qc_se_GD1.rds")
GD2 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_G2_C2/qc_se_GD2.rds")
GD3 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_G2_C2/qc_se_GD3.rds")

#Control
C1 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_G2_C2/qc_se_C1.rds")
C2 <- readRDS("~/Spatial_gene_expression/QC_spatial_level_G2_C2/qc_se_C2.rds")


# The name of the project
project_name <- "INTEGRATION_ALL_LogNorm_Harmony"

# 


# 
# ########
HT1@meta.data[["Condition"]] <- "HT"
HT2@meta.data[["Condition"]] <- "HT"
HT3@meta.data[["Condition"]] <- "HT"

GD1@meta.data[["Condition"]] <- "GD"
GD2@meta.data[["Condition"]] <- "GD"
GD3@meta.data[["Condition"]] <- "GD"


C1@meta.data[["Condition"]] <- "Controls"
C2@meta.data[["Condition"]] <- "Controls"

######### EXPERIMENT DAY
HT1@meta.data[["Experiment_day"]] <- "1"
HT2@meta.data[["Experiment_day"]] <- "1"
HT3@meta.data[["Experiment_day"]] <- "1"

GD1@meta.data[["Experiment_day"]] <- "1"
GD2@meta.data[["Experiment_day"]] <- "2"
GD3@meta.data[["Experiment_day"]] <- "2"


C1@meta.data[["Experiment_day"]] <- "2"
C2@meta.data[["Experiment_day"]] <- "2"

#Only infiltration

# Make a list with them
list_of_combinations <- list(HT1, HT2, HT3,
                             GD1, GD2, GD3,
                             C1, C2)

# Clustering resolutions in a vector
resolutions_to_clusters <- c(0.1, 0.3, 0.5, 0.75, 1)

# Where is the reactome?
Reactome <- read.gmt("c2.cp.reactome.v7.2.entrez.gmt")

# Cores to use
cores <- 15


# ================================================================================= #
##################################### DIRECTORY #####################################
# ================================================================================= #

dir.create(paste0("./Integration_", project_name))
wd_ana <- paste0("./Integration_", project_name)
setwd(wd_ana)


# ================================================================================= #
################################### EXTRA FUNCTIONS #################################
# ================================================================================= #

sumfun<-function(y){
  sum<-0
  for(i in 1:(length(y))){
    sum=sum+y[i]
  }
  print(sum)
}


# ================================================================================= #
######################## NORMALIZATION AND VARIABLE FEATURES ########################
# ================================================================================= #

list_of_combinations <- lapply(list_of_combinations, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, verbose = F)
})


var.features <- SelectIntegrationFeatures(object.list = list_of_combinations,
                                          nfeatures = 3000)

together <- merge(x = list_of_combinations[[1]],
                  y = list_of_combinations[2:length(list_of_combinations)],
                  merge.data=TRUE)

VariableFeatures(together) <- unique(var.features)


# ================================================================================= #
#################################### INTEGRATION ####################################
# ================================================================================= #

together<- ScaleData(together)
together <- RunPCA(together, verbose = FALSE)

together <- RunHarmony(together,features = var.features ,assay.use="Spatial",
                       group.by.vars = c("sample_id", "Experiment_day"))

ElbowPlot(together, reduction = "harmony")

# ================================================================================= #
############################## NEIGHBORS AND CLUSTERS ###############################
# ================================================================================= #


together <- RunUMAP(together, reduction = "harmony", dims = 1:13)
together <- FindNeighbors(together, reduction = "harmony", dims = 1:13)

together <- FindClusters(together, resolution = resolutions_to_clusters,
                         verbose = F)


# ================================================================================= #
################################# FIRST VIEW PLOTS ##################################
# ================================================================================= #

############## UMAP

#Create a vector with all the resolutions we have set
res_names <- paste("Spatial_snn_res", resolutions_to_clusters,
                   sep = ".")
#


plotPrinted<- DimPlot(object = together,
                      group.by = res_names,
                      label = TRUE,
                      label.size = 5,
                      reduction = "umap",
                      pt.size = 1.5)

#save
png(filename = paste0("./", project_name,"_","UMAP_diff_res.png"),
    width = 1000,
    height = 700)
print(plotPrinted)
dev.off()

############### PER SAMPLE

png(filename = paste0("./",project_name,"_persample_UMAP.png"),
    width = 1000,
    height = 700)
DimPlot(together, group.by = "sample_id",
        label.size = 10, pt.size = 1.5#,
        #cols = c("red", "blue", "gold3")
)
dev.off()

############### PER condition

png(filename = paste0("./",project_name,"_percondition_UMAP.png"),
    width = 1000,
    height = 700)
DimPlot(together, group.by = "Condition",
        label.size = 10, pt.size = 1.5#,
        #cols = c("red", "blue", "gold3")
)
dev.off()

############### CLUSTERS' TREE
png(filename = paste0("./",project_name,"_tree_clusters.png"),
    width = 1000,
    height = 700)
clustree(together)
dev.off()


############### CLUSTERS SPATIALLY REPRESENTED
n_to_w <- c("01", "03", "05", "075", "1")
pos <- 1
for (n_res in res_names){
  Idents(together)<- n_res
  plotPrinted <- SpatialDimPlot(together,
                                label = TRUE,
                                label.size = 5,
                                repel = T,
                                alpha = 0.9)
  
  png(filename = paste0("./", project_name,"_","clusters_res_",pos,".png"),
      width = 1000,
      height = 700)
  print(plotPrinted)
  dev.off()
  
  pos <- pos+1
  
}
remove(n_to_w, pos)
  


########### Removing iteratively a sample to evaluate integration at 1 res

img <- DimPlot(together, group.by = "Spatial_snn_res.0.5",
               label.size = 10, pt.size = 1.5)
png(filename = paste0("./",project_name,"_UMAP_0.5_",
                      ".png"),
    width = 1000,
    height = 700)
print(img)
dev.off()

Idents(together) <- "sample_id"
the_samples <- unique(together$sample_id)

for(i in c(1:length(the_samples))){
  sub_obj <- subset(together, idents = the_samples[-i])
  img <- DimPlot(sub_obj, group.by = "Spatial_snn_res.0.5",
                 label.size = 10, pt.size = 1.5)
  png(filename = paste0("./",project_name,"_UMAP_0.5_NO_", the_samples[i],
                        ".png"),
      width = 1000,
      height = 700)
  print(img)
  dev.off()
  
}

############# PCA
Idents(together) <- "Spatial_snn_res.0.5"
eigValues <- (together@reductions$pca@stdev)^2  ## EigenValues
varExplained <- eigValues / sum(eigValues)
plotPrinted <- DimPlot(together, label = T ,reduction = "pca",
                       label.size = 5)+
  xlab(paste0("PC1: ",round(varExplained[1]*100,1),"%"))+
  ylab(paste0("PC2: ",round(varExplained[2]*100,1),"%"))

#save
png(filename = paste0("./", project_name,"_","PCA_res.png"),
    width = 700,
    height = 700)
print(plotPrinted)
dev.off()

png(filename = paste0("./",project_name,"_percondition_pca.png"),
    width = 1000,
    height = 700)
DimPlot(together, group.by = "Condition",
        label.size = 10, pt.size = 1.5, reduction = "pca"#,
        #cols = c("red", "blue", "gold3")
)
dev.off()

############### PER condition

png(filename = paste0("./",project_name,"_pathologybased_UMAP.png"),
    width = 1000,
    height = 700)
DimPlot(together, group.by = "path_based",
        label.size = 10, pt.size = 1.5#,
        #cols = c("red", "blue", "gold3")
)
dev.off()


# ==================================================================== #
################################# MARKERS ##############################
# ==================================================================== #

dir.create("./Markers")

  
  for (res in res_names[3:5]){
    #Create directory for the resolution
    print("==================================================================")
    print(paste0("LOOP RESOLUTION: ", res))
    dir.create(paste0("./Markers/", res,"/"))
    to_save <- paste0("./Markers/", res,"/")
    #Which resolution?
    Idents(object = together) <- res
    #We obtain the markers
    markers <- FindAllMarkers(together)
    
    
    for(num in levels(markers$cluster)){
      #### Divide per cluster ####
      #### Markers
      markers_clust <- markers[which(markers$cluster==num),]
      gene_name <- mapIds(org.Hs.eg.db, keys=markers_clust$gene,
                          column="GENENAME", keytype="SYMBOL", multiVals="first")
      markers_clust["Gene_name"] <- gene_name
      #### GO
      sigOE <- dplyr::filter(as.data.frame(markers_clust), p_val_adj < 0.05)
      
      ego <- enrichGO(gene = as.character(rownames(sigOE)), 
                      #universe = allOE_genes,
                      keyType = "SYMBOL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = FALSE #TRUE when keytype is not SYMBOL
      )
      cluster_summary <- data.frame(ego)
      
      #### GSEA ####
      
      #### Foldchanges
      foldchanges <- markers_clust$avg_log2FC

      ## We have to transform Ensemble IDS to Entrez ID
      entrez_id = mapIds(org.Hs.eg.db, keys=markers_clust$gene,
                         column="ENTREZID", keytype="SYMBOL", multiVals="first")
      
      ## Name each fold change with the corresponding Entrez I
      names(foldchanges) <- entrez_id
      
      ## Sort fold changes in decreasing order
      foldchanges <- sort(foldchanges, decreasing = TRUE)
      
      #### Foldchanges end
      
      if (length(foldchanges) >50){
        gseaGO <- gseGO(geneList = foldchanges, 
                        OrgDb = org.Hs.eg.db, #to load the package
                        ont = 'BP',
                        nPerm = 1000, 
                        minGSSize = 20, 
                        pvalueCutoff = 0.05,
                        verbose = FALSE) 
        
        gseaGO_results <- gseaGO@result
        
        if ("core_enrichment" %in% colnames(gseaGO_results)){
          for (row in 1:nrow(gseaGO_results)){
            new_vect <- as.vector(unlist(strsplit(gseaGO_results[row,"core_enrichment"], "/")))
            symb_genes <- mapIds(org.Hs.eg.db, keys=new_vect,
                                 column="SYMBOL", keytype="ENTREZID", multiVals="first")
            gseaGO_results[row,"genes_in_path"] <- paste(symb_genes,collapse="/")
          }
          
        }
        
        
        #### KEGG ####
        gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                            organism = "hsa", # supported organisms listed below
                            nPerm = 1000, # default number permutations
                            #keyType = "ncbi-geneid",
                            minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                            pvalueCutoff = 0.05, # padj cutoff value
                            verbose = FALSE)
        
        
        gseaKEGG_results <- gseaKEGG@result
        
        if ("core_enrichment" %in% colnames(gseaKEGG_results)){
          for (row in 1:nrow(gseaKEGG_results)){
            new_vect <- as.vector(unlist(strsplit(gseaKEGG_results[row,"core_enrichment"], "/")))
            symb_genes <- mapIds(org.Hs.eg.db, keys=new_vect,
                                 column="SYMBOL", keytype="ENTREZID", multiVals="first")
            gseaKEGG_results[row,"genes_in_path"] <- paste(symb_genes,collapse="/")
          }
          
        }
        
        
        #### Reactome ####
        Reactome_sig <- GSEA(foldchanges, TERM2GENE=Reactome, verbose=FALSE)
        
        Reactome_sig_df <- data.frame(Reactome_sig)
        
        if ("core_enrichment" %in% colnames(Reactome_sig_df)){
          for (row in 1:nrow(Reactome_sig_df)){
            new_vect <- as.vector(unlist(strsplit(Reactome_sig_df[row,"core_enrichment"], "/")))
            symb_genes <- mapIds(org.Hs.eg.db, keys=new_vect,
                                 column="SYMBOL", keytype="ENTREZID", multiVals="first")
            Reactome_sig_df[row,"genes_in_path"] <- paste(symb_genes,collapse="/")
          }
          
        }
        
        #### Excel ####
        wb <- createWorkbook()
        #add cond
        addWorksheet(wb, sheetName="Differential_gene_expression")
        writeData(wb,sheet = "Differential_gene_expression", markers_clust)
        #add GO
        addWorksheet(wb, sheetName="GO")
        writeData(wb,sheet = "GO", cluster_summary)
        #add gseaGO
        addWorksheet(wb, sheetName="GSEA_GO")
        writeData(wb,sheet = "GSEA_GO", gseaGO_results)
        #add KEGG
        addWorksheet(wb, sheetName="KEGG")
        writeData(wb,sheet = "KEGG", gseaKEGG_results)
        #add reactome
        addWorksheet(wb, sheetName="Reactome")
        writeData(wb,sheet = "Reactome", Reactome_sig_df)
        #save
        saveWorkbook(wb, paste0(to_save ,res,"_cluster",num,".xlsx"))
        #### Plot ####
        plot <- SpatialFeaturePlot(object = together,
                                   features = markers_clust$gene[1:10],
                                   alpha = c(0.1, 1), ncol = 5)
        png(filename = paste0(to_save ,res,"_clust_",num, ".png"),
            width = 1000,height = 1000)
        print(plot)
        dev.off()
      }else{
        #### Excel ####
        wb <- createWorkbook()
        #add cond
        addWorksheet(wb, sheetName="Differential_gene_expression")
        writeData(wb,sheet = "Differential_gene_expression", markers_clust)
        #add GO
        addWorksheet(wb, sheetName="GO")
        writeData(wb,sheet = "GO", cluster_summary)
        #save
        saveWorkbook(wb, paste0(to_save ,res,"_cluster",num,".xlsx"))
        #### Plot ####
        plot <- SpatialFeaturePlot(object = together,
                                   features = markers_clust$gene[1:10],
                                   alpha = c(0.1, 1), ncol = 5)
        png(filename = paste0(to_save ,res,"_clust_",num, ".png"),
            width = 1000,height = 1000)
        print(plot)
        dev.off()
      }

      
    }
    
  }
  

  

# ================================================================================= #
######################################## SAVE #######################################
# ================================================================================= #

save.image(paste0("./",project_name, "_R_environment.RData"))

