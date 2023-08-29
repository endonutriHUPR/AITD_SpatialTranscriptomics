# Re-clustering

Idents(together) <- "path_based"
cell_population <- subset(x = together, idents = c("cell_population"))
cell_population_list <- SplitObject(cell_population, split.by = "sample_id")
cell_population_list <- lapply(X = cell_population_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
int.features <- SelectIntegrationFeatures(object.list = cell_population_list,
                                          nfeatures = 2000)


cell_population_reclust <- merge(x = cell_population_list[[1]],
                  y = thyro_list[2:length(thyro_list)],
                  merge.data=TRUE)
remove(cell_population_list,cell_population )

VariableFeatures(cell_population_reclust) <- unique(int.features)

cell_population_reclust<- ScaleData(cell_population_reclust)
cell_population_reclust <- RunPCA(cell_population_reclust, verbose = FALSE, npcs = 30)

## Re-integration
cell_population_reclust <- RunHarmony(cell_population_reclust,dims.use = 1:30,assay.use = "Spatial",
                       group.by.vars = c("sample_id",  "Experiment_day"))
                       
                       
ElbowPlot(cell_population_reclust, reduction = "harmony")


cell_population_reclust <- RunUMAP(cell_population_reclust, reduction = "harmony", dims = 1:12)


cell_population_reclust <- FindNeighbors(cell_population_reclust, reduction = "harmony", dims = 1:12)


cell_population_reclust <- FindClusters(cell_population_reclust, resolution = c(0.1, 0.3, 0.5, 0.75, 1),
                         verbose = F)
Idents(cell_population_reclust) <- "sample_id"
DimPlot(object = cell_population_reclust,
                      label = FALSE, label.size = 5,pt.size = 0.75,
                      reduction = "umap")&
                     paletteer::scale_colour_paletteer_d("tidyquant::tq_light")
                     
                     
Idents(cell_population_reclust) <- "Condition"
DimPlot(object = cell_population_reclust,
                      label = FALSE, label.size = 5,pt.size = 0.75,
                      reduction = "umap")&
                     paletteer::scale_colour_paletteer_d("tidyquant::tq_light")
                     
                     
                     
Idents(cell_population_reclust) <- "Spatial_snn_res.0.1"
pt1 <- DimPlot(object = cell_population_reclust,
                      #label.size = 5,
        pt.size = 1.25,
                      reduction = "umap",

                    repel = T
                    )&
                     paletteer::scale_colour_paletteer_d("tidyquant::tq_light")

Idents(cell_population_reclust) <- "Spatial_snn_res.0.3"

pt2<- DimPlot(object = cell_population_reclust,
                      #label.size = 5,
        pt.size = 1.25,
                      reduction = "umap",

                    repel = T
                    )&
                     paletteer::scale_colour_paletteer_d("tidyquant::tq_light")

Idents(cell_population_reclust) <- "Spatial_snn_res.0.5"

pt3<-DimPlot(object = cell_population_reclust,
                      #label.size = 5,
        pt.size = 1.25,
                      reduction = "umap",

                    repel = T
                    )&
                     paletteer::scale_colour_paletteer_d("tidyquant::tq_light")
cowplot::plot_grid(pt1,pt2,pt3, labels = c("Res 0.1", "Res 0.3", "Res 0.5"), label_size = 20)


Idents(cell_population_reclust) <- "Spatial_snn_res.0.5"
cell_population_reclust_markers <- FindAllMarkers(cell_population_reclust)



