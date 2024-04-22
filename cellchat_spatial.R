library(CellChat)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 
library(future)
library(future.apply)


CellChatDB <- CellChatDB.human


#HT
data.input <- Seurat::GetAssayData(HT, slot = "data", assay = "Spatial")
meta = data.frame(labels = Idents(HT), samples = HT@meta.data$sample_id) # manually create a dataframe consisting of the cell labels
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(HT)))
meta$samples <- factor(meta$samples, levels = unique(HT@meta.data$sample_id))
meta$labels = droplevels(meta$labels, exclude = NA)

spatial.locs <- Seurat::GetTissueCoordinates(HT, scale = NULL, cols = c("imagerow", "imagecol"), image = "HT1") 
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(HT, scale = NULL, cols = c("imagerow", "imagecol"), image = "HT2") )
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(HT, scale = NULL, cols = c("imagerow", "imagecol"), image = "HT3") )

spot.size = 65 # the theoretical spot size (um) in 10X Visium

#list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres
scalefactors1 = jsonlite::fromJSON(txt = file.path("/HT1", 'scalefactors_json.json'))
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors1$spot_diameter_fullres
                              )

scalefactors2 = jsonlite::fromJSON(txt = file.path("/HT2", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors2$spot_diameter_fullres
                              )

scalefactors3 = jsonlite::fromJSON(txt = file.path("/HT3", 'scalefactors_json.json'))
conversion.factor3 = spot.size/scalefactors3$spot_diameter_fullres
spatial.factors3 = data.frame(ratio = conversion.factor3, tol = spot.size/2#,
                              #spot.diameter = spot.size, spot=scalefactors3$spot_diameter_fullres
                              )

spatial.factors <- rbind(spatial.factors1, spatial.factors2,
                         spatial.factors3)

rownames(spatial.factors) <- unique(HT@meta.data$sample_id)

remove(scalefactors1, scalefactors2, scalefactors3)


cellchat_HT_ALL <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs,spatial.factors = spatial.factors)
                           
cellchat_HT_ALL@DB <- CellChatDB
cellchat_HT_ALL <- subsetData(cellchat_HT_ALL) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) 
options(future.rng.onMisuse="ignore")

cellchat_HT_ALL <- identifyOverExpressedGenes(cellchat_HT_ALL)
cellchat_HT_ALL <- identifyOverExpressedInteractions(cellchat_HT_ALL)


cellchat_HT_ALL<- computeCommunProb(cellchat_HT_ALL, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250,
                              scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100
                              )

cellchat_HT_ALL <- computeCommunProbPathway(cellchat_HT_ALL)

cellchat_HT_ALL <- filterCommunication(cellchat_HT_ALL, min.cells = 10)


cellchat_HT_ALL <- computeCommunProbPathway(cellchat_HT_ALL)
cellchat_HT_ALL <- aggregateNet(cellchat_HT_ALL)


# GD

data.input <- Seurat::GetAssayData(GD, slot = "data", assay = "Spatial")
meta = data.frame(labels = Idents(GD), samples = GD@meta.data$sample_id) # manually create a dataframe consisting of the cell labels
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(HT)))
meta$samples <- factor(meta$samples, levels = unique(GD@meta.data$sample_id))
meta$labels = droplevels(meta$labels, exclude = NA)

spatial.locs <- Seurat::GetTissueCoordinates(GD, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD1") 
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(GD, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD2") )
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(GD, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD3") )

spot.size = 65 # the theoretical spot size (um) in 10X Visium

#list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres
scalefactors1 = jsonlite::fromJSON(txt = file.path("/GD1", 'scalefactors_json.json'))
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors1$spot_diameter_fullres
                              )

scalefactors2 = jsonlite::fromJSON(txt = file.path("/GD2", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors2$spot_diameter_fullres
                              )

scalefactors3 = jsonlite::fromJSON(txt = file.path("/GD3", 'scalefactors_json.json'))
conversion.factor3 = spot.size/scalefactors3$spot_diameter_fullres
spatial.factors3 = data.frame(ratio = conversion.factor3, tol = spot.size/2#,
                              #spot.diameter = spot.size, spot=scalefactors3$spot_diameter_fullres
                              )

spatial.factors <- rbind(spatial.factors1, spatial.factors2,
                         spatial.factors3)

rownames(spatial.factors) <- unique(GD@meta.data$sample_id)

remove(scalefactors1, scalefactors2, scalefactors3)


cellchat_GD_ALL <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs,spatial.factors = spatial.factors)
                           
cellchat_GD_ALL@DB <- CellChatDB
cellchat_GD_ALL <- subsetData(cellchat_GD_ALL) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) 
options(future.rng.onMisuse="ignore")

cellchat_GD_ALL <- identifyOverExpressedGenes(cellchat_GD_ALL)
cellchat_GD_ALL <- identifyOverExpressedInteractions(cellchat_GD_ALL)


cellchat_GD_ALL<- computeCommunProb(cellchat_GD_ALL, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250,
                              scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100
                              )

cellchat_GD_ALL <- computeCommunProbPathway(cellchat_GD_ALL)

cellchat_GD_ALL <- filterCommunication(cellchat_GD_ALL, min.cells = 10)


cellchat_GD_ALL <- computeCommunProbPathway(cellchat_GD_ALL)
cellchat_GD_ALL <- aggregateNet(cellchat_GD_ALL)

# Controls

data.input <- Seurat::GetAssayData(CONTROLS, slot = "data", assay = "Spatial")
meta = data.frame(labels = Idents(CONTROLS), samples = CONTROLS@meta.data$sample_id) # manually create a dataframe consisting of the cell labels
rownames(meta) <- colnames(data.input)
meta$labels <- factor(meta$labels, levels = levels(Idents(CONTROLS)))
meta$samples <- factor(meta$samples, levels = unique(CONTROLS@meta.data$sample_id))
meta$labels = droplevels(meta$labels, exclude = NA)

spatial.locs <- Seurat::GetTissueCoordinates(CONTROLS, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD1") 
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(CONTROLS, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD2") )
spatial.locs <- rbind(spatial.locs, Seurat::GetTissueCoordinates(CONTROLS, scale = NULL, cols = c("imagerow", "imagecol"), image = "GD3") )

spot.size = 65 # the theoretical spot size (um) in 10X Visium

#list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres
scalefactors1 = jsonlite::fromJSON(txt = file.path("/CONTROLS1", 'scalefactors_json.json'))
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors1$spot_diameter_fullres
                              )

scalefactors2 = jsonlite::fromJSON(txt = file.path("/CONTROLS2", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2#,
                              #spot.diameter = spot.size,spot=scalefactors2$spot_diameter_fullres
                              )

scalefactors3 = jsonlite::fromJSON(txt = file.path("/CONTROLS3", 'scalefactors_json.json'))
conversion.factor3 = spot.size/scalefactors3$spot_diameter_fullres
spatial.factors3 = data.frame(ratio = conversion.factor3, tol = spot.size/2#,
                              #spot.diameter = spot.size, spot=scalefactors3$spot_diameter_fullres
                              )

spatial.factors <- rbind(spatial.factors1, spatial.factors2,
                         spatial.factors3)

rownames(spatial.factors) <- unique(CONTROLS@meta.data$sample_id)

remove(scalefactors1, scalefactors2, scalefactors3)


cellchat_CONTROLS_ALL <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs,spatial.factors = spatial.factors)
                           
cellchat_CONTROLS_ALL@DB <- CellChatDB
cellchat_CONTROLS_ALL <- subsetData(cellchat_CONTROLS_ALL) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) 
options(future.rng.onMisuse="ignore")

cellchat_CONTROLS_ALL <- identifyOverExpressedGenes(cellchat_CONTROLS_ALL)
cellchat_CONTROLS_ALL <- identifyOverExpressedInteractions(cellchat_CONTROLS_ALL)


cellchat_CONTROLS_ALL<- computeCommunProb(cellchat_CONTROLS_ALL, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250,
                              scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100
                              )

cellchat_CONTROLS_ALL <- computeCommunProbPathway(cellchat_CONTROLS_ALL)

cellchat_CONTROLS_ALL <- filterCommunication(cellchat_CONTROLS_ALL, min.cells = 10)


cellchat_CONTROLS_ALL <- computeCommunProbPathway(cellchat_CONTROLS_ALL)
cellchat_CONTROLS_ALL <- aggregateNet(cellchat_CONTROLS_ALL)


