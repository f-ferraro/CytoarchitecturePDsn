#####################################################################
#################### Preprocessing of GSE140231 #####################
#####################################################################

# remotes::install_version("Seurat", version = "2.3.0")

library(dplyr)
library(Seurat)
library(patchwork)
library(GEOquery)

meta <-  read.csv("../2resources/Nigra/GSE140231_meta.csv")

GSM4157078 <- Read10X(data.dir = "../2resources/Nigra/GSM4157078")
GSM4157076 <- Read10X(data.dir = "../2resources/Nigra/GSM4157076")
GSM4157074 <- Read10X(data.dir = "../2resources/Nigra/GSM4157074")
GSM4157072 <- Read10X(data.dir = "../2resources/Nigra/GSM4157072")
GSM4157070 <- Read10X(data.dir = "../2resources/Nigra/GSM4157070")
GSM4157069 <- Read10X(data.dir = "../2resources/Nigra/GSM4157069")
GSM4157068 <- Read10X(data.dir = "../2resources/Nigra/GSM4157068")


list <- list(GSM4157078, GSM4157076, GSM4157074, GSM4157072,
             GSM4157070, GSM4157069, GSM4157068)
names(list) <- c("GSM4157078", "GSM4157076", "GSM4157074", "GSM4157072",
                 "GSM4157070", "GSM4157069", "GSM4157068")
for (proband in names(list)){
  obj<- list[[proband]]
  obj <- CreateSeuratObject(obj,min.cells = 3, min.features = 500)
  mito.genes <- grep("^MT-", rownames(obj@data), value = T)
  percent.mito <- colSums(expm1(obj@data[mito.genes, ]))/colSums(expm1(obj@data))
  RPS.genes <- grep(pattern = "^RPS", x = rownames(x = obj@data), value = TRUE)
  percent.RPS <- Matrix::colSums(obj@raw.data[RPS.genes, ])/Matrix::colSums( obj@raw.data)
  obj <- AddMetaData(object = obj, metadata = percent.RPS, col.name = "percent.RPS")
  obj <- AddMetaData(obj, percent.mito, "percent.mito")
  # VlnPlot(obj, c("nGene", "nUMI", "percent.mito", "percent.RPS"), nCol = 2)
  # par(mfrow = c(1, 2))
  # GenePlot(obj, "nUMI", "percent.mito")
  # GenePlot(obj, "nUMI", "nGene")
  # GenePlot(obj, "nUMI", "percent.RPS")
  obj <- SubsetData(obj, subset.name = "nGene", accept.high = 5000)
  obj <- SubsetData(obj, subset.name = "percent.mito", accept.high = 0.05)
  obj <- SubsetData(obj, subset.name = "percent.RPS", accept.high = 0.05)
  obj <- NormalizeData(object = obj, normalization.method = "LogNormalize",
                       scale.factor = 10000)
  obj <- ScaleData(object = obj, vars.to.regress = c("nUMI",
                                                     "percent.mito",
                                                     "percent.RPS"))


  list[[proband]] <- obj
}

# save.image("../2resources/GSMquality.RData")
load("../2resources/GSMquality.RData")

rm(list=setdiff(ls(), "list"))

list2env(list, globalenv())

SNC <- MergeSeurat(object1 = GSM4157078, object2 = GSM4157076, add.cell.id1 = "GSM4157078",
                   add.cell.id2 = "GSM4157076", project = "SN", do.normalize = FALSE)

SNC3 <- MergeSeurat(object1 = GSM4157070, object2 = GSM4157072, add.cell.id1 = "GSM4157070",
                    add.cell.id2 = "GSM4157072", project = "SN", do.normalize = FALSE)

SNC4 <- MergeSeurat(object1 = GSM4157068, object2 = GSM4157069, add.cell.id1 = "GSM4157068",
                    add.cell.id2 = "GSM4157069", project = "SN", do.normalize = FALSE)


SNC5 <- MergeSeurat(object1 = SNC, object2 = GSM4157074, add.cell.id1 = "",
                    add.cell.id2 = "GSM4157074", project = "SN", do.normalize = FALSE)
SNC6 <- MergeSeurat(object1 = SNC3, object2 = SNC4, add.cell.id1 = "",
                    add.cell.id2 = "", project = "SN", do.normalize = FALSE)

SN1 <- MergeSeurat(object1 = SNC5, object2 = SNC6, add.cell.id1 = "",
                   add.cell.id2 = " ", project = "SN", do.normalize = FALSE)

rm(list=setdiff(ls(), "SN1"))

SNC <- SN1
rm(list=setdiff(ls(), "SNC"))
SNC <- NormalizeData(object = SNC)
SNC <- FindVariableGenes(object = SNC, mean.function = ExpMean, dispersion.function = LogVMR,
                         x.low.cutoff = 0.02, x.high.cutoff = 3.5, y.cutoff = 0.5)

length(x = SNC@var.genes)

SNC <- ScaleData(SNC, genes.use = SNC@var.genes)

SNC <- RunPCA(object = SNC, pc.genes = SNC@var.genes, pcs.compute = 40,
              do.print = FALSE, pcs.print = 1:5, genes.print = 5)

# SNC <- JackStraw(object = SNC, num.replicate = 100, display.progress = FALSE)
# JackStrawPlot(object = SNC, PCs = 1:30)
PCElbowPlot(object = SNC, num.pc = 40)

#####################################################################
################ Markers identification GSE140231 ###################
#####################################################################
SNC <- FindClusters(object = SNC, k.param=30, reduction.type = "pca", dims.use = 1:25,
                    resolution = 0.4, print.output = 0, save.SNN = TRUE)

#Plot
SNC <- RunTSNE(object = SNC, dims.use = 1:25, do.fast = TRUE)
TSNEPlot(object = SNC)

saveRDS(SNC, file = "../2resources/Tsne_SNC.rds")
# FeaturePlot(object = SNC, features.plot = c("VCAN"), cols.use = c("grey", "blue"),
# reduction.use = "tsne")

# find markers for every cluster compared to all remaining cells
SNC <- readRDS("~/Scrivania/github/2resources/Tsne_SNC.rds")
SNC.markers <- FindAllMarkers(object = SNC, only.pos = TRUE, min.pct = 0.30, 
                              thresh.use = 0.5)
# Markers from literature for cluster annotation
goi <- c("GFAP", "GINS3", "MOBP", "GINS3", "TH", "SLC6A3",  "GABRA1", "GABRA2", 
         "GAD1", "GAD2", "RGS5", "MOG", "OLR1", "CSF1R", "PPM1G", "LGALS1", "PALM2","VCAN")
goi <- unique(as.character(goi))
DotPlot(SNC,genes.plot = (goi))

# Look at the heatmap, choose the marker and rename 
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
new.cluster.ids <- c("ODC", #0
                     "ODC", #1
                     "ODC", #2 
                     "Astrocyte", #3
                     "Microglia", #4
                     "OPC", #5  
                     "Neuron", #6DaNs
                     "Endothelial", #7 
                     "Neuron") #8 GABA  
names(new.cluster.ids) <- levels(SNC)
SNC@ident <- plyr::mapvalues(x = SNC@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = SNC, do.label = TRUE, pt.size = 0.5)

DoHeatmap(SNC,genes.use = goi,slim.col.label = TRUE,remove.key = TRUE)


#######################################################################

SNC.markers$cluster <- gsub(0, "ODC", SNC.markers$cluster)
SNC.markers$cluster <- gsub(1, "ODC", SNC.markers$cluster)
SNC.markers$cluster <- gsub(2, "ODC", SNC.markers$cluster)
SNC.markers$cluster <- gsub(3, "Astrocyte", SNC.markers$cluster)
SNC.markers$cluster <- gsub(4, "Microglia", SNC.markers$cluster)
SNC.markers$cluster <- gsub(5, "OPC", SNC.markers$cluster)
SNC.markers$cluster <- gsub(6, "DaNs", SNC.markers$cluster)
SNC.markers$cluster <- gsub(7, "Endothelial", SNC.markers$cluster)
SNC.markers$cluster <- gsub(8, "GABA", SNC.markers$cluster)


names(new.cluster.ids) <- levels(SNC)


write.csv(SNC.markers, "../2resources/SNC.markers")