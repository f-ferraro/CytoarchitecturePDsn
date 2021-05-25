##############################################################################################
############################### DEGs_EWCE in GSE140231 #######################################
##############################################################################################

library("EWCE")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("GEOquery")



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

SNC <- FindClusters(object = SNC, k.param=30, reduction.type = "pca", dims.use = 1:25,
                    resolution = 0.4, print.output = 0, save.SNN = TRUE)

SNC <- RunTSNE(object = SNC, dims.use = 1:25, do.fast = TRUE)
TSNEPlot(object = SNC)
SNC.markers <- FindAllMarkers(object = SNC, only.pos = TRUE, min.pct = 0.30,
                              thresh.use = 0.5)
goi <- c("GFAP", "GINS3", "MOBP", "GINS3", "TH", "SLC6A3",  "GABRA1", "GABRA2",
         "GAD1", "GAD2", "RGS5", "MOG", "OLR1", "CSF1R", "PPM1G", "LGALS1", "PALM2","VCAN")
goi <- unique(as.character(goi))
DotPlot(SNC,genes.plot = (goi))

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

l1 <- as.list(SNC@ident)
l1 <- lapply(l1, as.character)
l1 <- as.character(SNC@ident)
l1 <- list(l1, l1)
exp <- as.data.frame(SNC@scale.data)
gcd <- generate.celltype.data(exp,
                              l1,
                              "SNCA",
                              )

mydf <- read.delim2("../2resources/averageexpression.tsv", row.names = 1)
# mydf <- log2(mydf)
# mydf[mydf=="-Inf"] <- 0
l1 <- as.character(names(mydf))
l1 <- list(l1, l1)
gcd <- generate.celltype.data(mydf,
                              l1,
                              "SN")
load("../2resources/CellTypeData_SN.rda")
corrected <- read.csv2("~/Scrivania/eNeuro/3results/correctedDEGs.csv", row.names=1)
corrected$HGNC.symbol <- corrected$Gene
corrected <- corrected[corrected$pBH<0.05 &abs(corrected$estimate)>log(1.2),]
corrected <- corrected[corrected$HGNC.symbol %in% rownames(mydf),]
corrected$Gene <- as.character(corrected$Gene)
rm(list=setdiff(ls(), c("dnDEG", "upDEG", "ctd","corrected")))

##############################################################
# EWCE plot DEGs
full_res_up = EWCE::bootstrap.enrichment.test(sct_data = ctd, 
                                              hits = corrected[corrected$estimate>log(1.2) &
                                                                 corrected$pBH<0.05,]$Gene, 
                                              bg = rownames(ctd[[1]][["mean_exp"]]),
                                              reps = 10000, 
                                              annotLevel = 1, 
                                              geneSizeControl = FALSE, 
                                              genelistSpecies = "human", 
                                              sctSpecies = "human")

full_res_down = EWCE::bootstrap.enrichment.test(sct_data = ctd, 
                                                hits = corrected[corrected$estimate<log(1.2) &
                                                                   corrected$pBH<0.05,]$Gene, 
                                                bg = rownames(ctd[[1]][["mean_exp"]]),
                                                reps = 10000, 
                                                annotLevel = 1, 
                                                geneSizeControl = FALSE, 
                                                genelistSpecies = "human", 
                                                sctSpecies = "human")

joint_results = rbind(cbind(full_res_up$results, Direction = "Up"), 
                      cbind(full_res_down$results, Direction = "Down"))
joint_results

joint_results$pBH <- p.adjust(joint_results$p, method = "BH")
joint_results$CellType <- gsub("ODC", "ODCs", joint_results$CellType)
joint_results$CellType <- gsub("Astrocyte", "Astrocytes", joint_results$CellType)
joint_results$CellType <- gsub("Endothelial", "Endothelial Cells", joint_results$CellType)
joint_results$CellType <- gsub("OPC", "OPCs", joint_results$CellType)
joint_results$CellType <- gsub("Neuron", "Neurons", joint_results$CellType)

joint_results$CellType <- gsub("Endothelial Cells", "Endothelial \nCells", joint_results$CellType)
tiff("../4plots/Figure3E.tiff", width = 9.5, height = 7, units = "in", res = 300)
ggplot(joint_results, aes(x=CellType, y=abs(sd_from_mean), fill=Direction, fill=Direction)) +
  geom_bar(stat='identity', position='dodge', alpha=0.8)+
  xlab("")+ylab("Std.Devs. from the mean")+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(labels = c("Up", "Down"),
                    values=c('#bb0021e5','#008b45e5')) +
  theme(legend.position = c(0, 1),text = element_text(size=20),
        legend.justification = c(-0.07,1),
        legend.title = element_text(color = "blue", size = 0))
dev.off()

joint_results

library(openxlsx)

write.xlsx(joint_results, file = "../3results/EDT3-1.xlsx")
#############################################################
# Plot heatmap preferential expression
mydf <- read.delim2("../2resources/averageexpression.tsv", row.names = 1)
corrected <- read.csv2("~/Scrivania/eNeuro/3results/correctedDEGs.csv", row.names=1)
corrected <- corrected[corrected$pBH<0.05 &abs(corrected$estimate)>log(1.2),]
corrected$Gene <- as.character(corrected$Gene)
corrected$HGNC.symbol <- corrected$Gene
corrected <- corrected[corrected$HGNC.symbol %in% rownames(mydf),]

rm(list=setdiff(ls(), c("mydf", "corrected")))
mydf <- mydf[rownames(mydf) %in% corrected$Gene,]
mydf <- t(scale(t(mydf)))

b <-corrected[, c("Gene", "estimate")]
rownames(b) <- b$Gene
b$Gene <- NULL
b$DEGs <- ifelse(b$estimate>0, "upDEGs", "dnDEGs")
b$estimate <- NULL

ann = list(DEGs=c(upDEGs="orangered3", dnDEGs="palegreen3"))


library(pheatmap)
library(RColorBrewer)
tiff("../4plots/SF3_preferentialexpression.tiff", height = 7, width = 5, res = 300, units = "in")
pheatmap(mydf,
         annotation_row=b,
         color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                   "RdBu")))(150),
         annotation_colors = ann, 
         show_rownames = F,
         treeheight_row = 0, treeheight_col = 0)
dev.off()


##############################################################################################
########################### GSEA LMMs-ranked gene expression  #################################
##############################################################################################
load("../3results/Original.RData")
mixed2 <- mixed
mixed2 <- bind_rows(mixed2, .id="Gene")
singfilt <- bind_rows(singular, .id="Gene")
mixed2 <- mixed2[!(mixed2$Gene %in% as.data.frame(t(singfilt))$V1), ]
mixed2 <- split(mixed2, as.factor(mixed2$term))
mixed2 <- mixed2$Status1
mixed2$pBH <- p.adjust(mixed2$p.value, method = "BH")
original <- mixed2
rm(list=setdiff(ls(), c("corrected","original", "resCB", "resCB2")))

load("../3results/Corrected.RData")
mixed2 <- mixed
mixed2 <- bind_rows(mixed2, .id="Gene")
singfilt <- bind_rows(singular, .id="Gene")
mixed2 <- mixed2[!(mixed2$Gene %in% as.data.frame(t(singfilt))$V1), ]
mixed2 <- split(mixed2, as.factor(mixed2$term))
mixed2 <- mixed2$Status1
mixed2$pBH <- p.adjust(mixed2$p.value, method = "BH")
corrected <- mixed2
rm(list=setdiff(ls(), c("corrected","original", "resCB", "resCB2")))

corrected = corrected[corrected$Gene %in% original$Gene, ]
original = original[original$Gene %in% corrected$Gene, ]

library("fgsea")
library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
canonical <- gmtPathways('../2resources/canonical_all.gmt')
GOBP <- gmtPathways('../2resources/GOBP.gmt')
GOMF <- gmtPathways('../2resources/GOMF.gmt')
GOCC <- gmtPathways('../2resources/GOCC.gmt')

annot <- Hmisc::llist(canonical,
                      GOBP, GOMF, GOCC)
rm(hallmark, curatedpath, canonical, GOBP, GOMF, GOCC, TF, mirna, position)

original$Score <- sign(original$estimate)*(-log10(original$pBH)) #(-log(original$pBH))
original <- dplyr::select(original, "Gene", "estimate", "pBH", "Score")
original <- original[order(original$Score, decreasing = T),]

entrez <- getBM(attributes= c("hgnc_symbol", "entrezgene_id"),
                values = original$Gene,
                mart = ensembl)
names(entrez) <- c("Gene", "ENTREZ" )
original <- merge(original, entrez, by="Gene")
original <- na.omit(original)
original <- original[order(-original$Score), ]
ranks_o <- original$Score
names(ranks_o) <- original$ENTREZ

#Corrected
corrected$Score <- sign(corrected$estimate)*(-log10(corrected$pBH))
corrected <- dplyr::select(corrected, "Gene", "estimate", "pBH", "Score")
corrected <- corrected[order(corrected$Score, decreasing = T),]

entrez <- getBM(attributes= c("hgnc_symbol", "entrezgene_id"),
                values = corrected$Gene,
                mart = ensembl)
names(entrez) <- c("Gene", "ENTREZ" )
corrected <- merge(corrected, entrez, by="Gene")
corrected <- na.omit(corrected)
corrected <- corrected[order(-corrected$Score), ]
ranks_c <- corrected$Score
names(ranks_c) <- corrected$ENTREZ

rm(list=setdiff(ls(), c("annot", "ranks_c", "ranks_o")))
gc()

corrected <- original <- list()
for (i in names(annot)){
  fgseaRes_c <- fgsea(annot[[i]], ranks_c, minSize = 25, maxSize=500, nperm=100000)
  fgseaRes_c$padj <- p.adjust(fgseaRes_c$pval, method = "BH")
  write.csv2(as.matrix(fgseaRes_c), paste0(i, "_corrected.csv"))
  
  
  fgseaRes_o <- fgsea(annot[[i]], ranks_o,  minSize = 25, maxSize=500, nperm=100000)
  fgseaRes_o$padj <- p.adjust(fgseaRes_o$pval, method = "BH")
  write.csv2(as.matrix(fgseaRes_o), paste0(i, "_original.csv"))
  
  corrected[[i]] <- fgseaRes_c
  original[[i]] <- fgseaRes_o
  rm(fgseaRes_o, fgseaRes_c)
  gc()}

save.image("../3results/functional.RData")
# load("../3results/functional.RData")

#############################################################################
# Change of significance 
GOBPc <- corrected$GOBP
GOCCc <- corrected$GOCC
GOMFc <- corrected$GOMF
GOBPo <- original$GOBP
GOCCo <- original$GOCC
GOMFo <- original$GOMF

db <- names(corrected)

for (i in names(corrected)) {
  names(corrected[[i]]) <-  paste0("corrected_",names(corrected[[i]]))
  names(original[[i]]) <- paste0("original_", names(original[[i]]))}

for (i in names(corrected)) {
  original[[i]]$path <- original[[i]]$original_pathway
  corrected[[i]]$path <- corrected[[i]]$corrected_pathway
}

compar <- list()
for (i in names(corrected)) {
  compar[[i]] <- merge(corrected[[i]], original[[i]], by = "path", all = TRUE)
  # compar[[i]]$corrected_padj <- -log10(compar[[i]]$corrected_padj)
  # compar[[i]]$original_padj <- -log10(compar[[i]]$original_padj)
  # compar[[i]] <- compar[[i]][compar[[i]]$corrected_padj>(-log10(0.05))|
  #                              compar[[i]]$original_padj>(-log10(0.05)), ]
  compar[[i]] <- as.data.frame(compar[[i]])
  # compar[[i]]$Delta <- compar[[i]]$original_padj - compar[[i]]$corrected_padj
  # compar[[i]]$Delta <- compar[[i]]$Delta/max(compar[[i]]$Delta, na.rm = T)
  # compar[[i]]$Nes <- ifelse(sign(compar[[i]]$corrected_NES)==sign(compar[[i]]$original_NES), "T", "F")
  # compar[[i]]$DNes <- (compar[[i]]$original_NES/compar[[i]]$corrected_NES)^2
  # compar[[i]]$pt <- compar[[i]]$path
  # compar[[i]]$pt <- gsub("_", " ", compar[[i]]$pt)
}

library(openxlsx)

for (i in names(compar)) {
  # compar[[i]]$corrected_padj <- 10^(-(compar[[i]]$corrected_padj))
  # compar[[i]]$original_padj <- 10^(-(compar[[i]]$original_padj))
  compar[[i]] <-  as.data.frame(compar[[i]])
  compar[[i]] <- dplyr::select(compar[[i]], c("path", "corrected_pval", "corrected_padj", "corrected_ES", "corrected_NES", "corrected_nMoreExtreme", "corrected_size",
                                              "original_pval", "original_padj", "original_ES", "original_NES",
                                              "original_nMoreExtreme", "original_size"))
}

write.xlsx(compar, file = "../3results/ST2_GSEA_comparison.xlsx")

rm(list = ls())
##############################################################################
# Plots
x <- read_excel("../3results/ST2_GSEA_comparison.xlsx", 
                sheet = "canonical")
names(x) <- gsub("CellAware", "corrected", names(x))
names(x) <- gsub("CellUnaware", "original", names(x))

library("ggplot2")
x <- x[order(x$corrected_padj),]
x <- x[1:10, ]
x$path <- gsub("_", " ", x$path)
x$path <- gsub("REACTOME", "REACTOME - ", x$path)
x$path <- gsub("KEGG", "KEGG - ", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$corrected_padj)])
tiff("../4plots/canonicalGSEA.tiff", units = "cm", width = 10, height = 12, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(corrected_padj)), data=x)+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(original_padj)), data=x)+
  geom_point( aes(x=path, y=-log10(original_padj), colour =original_NES), shape=16, size=4)+
  geom_point( aes(x=path, y=-log10(corrected_padj),  colour =corrected_NES), shape=18, size=4)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  theme(legend.position = "bottom",  legend.text = element_text(size = 6), 
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed")+
  labs(colour = "NES") + 
  scale_color_gradient2(low = "#008b45e5", high = "#bb0021e5",  limit = c(-2.5,2.5))+
  ggtitle("Canonical data set")
dev.off()

############################################################################################
library(readxl)
x <- read_excel("~/Desktop/eNeuro/3results/ST2_GSEA_comparison.xlsx", 
                sheet = "GOBP")
names(x) <- gsub("CellAware", "corrected", names(x))
names(x) <- gsub("CellUnaware", "original", names(x))
x <- x[order(x$corrected_padj),]
x <- x[1:10, ]
x$path <- gsub("_", " ", x$path)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$corrected_padj)])
tiff("../4plots/GOBP_GSEA.tiff", units = "cm", width = 10, height = 12, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(corrected_padj)), data=x)+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(original_padj)), data=x)+
  geom_point( aes(x=path, y=-log10(original_padj), colour =original_NES), shape=16, size=4)+
  geom_point( aes(x=path, y=-log10(corrected_padj),  colour =corrected_NES), shape=18, size=4)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  ggtitle("GOBP")+
  theme(legend.position = "bottom",  legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed")+ labs(colour = "NES")  + 
  scale_color_gradient2(low = "#008b45e5", high = "#bb0021e5",  limit = c(-2.5,2.5))

dev.off()
############################################################################################
library(readxl)
x <- read_excel("~/Desktop/eNeuro/3results/ST2_GSEA_comparison.xlsx", 
                sheet = "GOMF")
names(x) <- gsub("CellAware", "corrected", names(x))
names(x) <- gsub("CellUnaware", "original", names(x))
x <- x[order(x$corrected_padj),]
x <- x[1:10, ]
x$path <- gsub("_", " ", x$path)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$corrected_padj)])
tiff("../4plots/GOMF_GSEA.tiff", units = "cm", width = 10, height = 12, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(corrected_padj)), data=x)+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(original_padj)), data=x)+
  geom_point( aes(x=path, y=-log10(original_padj), colour =original_NES), shape=16, size=4)+
  geom_point( aes(x=path, y=-log10(corrected_padj),  colour =corrected_NES), shape=18, size=4)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic()  +
  ggtitle("GOMF")+
  theme(legend.position = "bottom",  legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed")+labs(colour = "NES") + 
  scale_color_gradient2(low = "#008b45e5", high = "#bb0021e5",  limit = c(-2.5,2.5)) 

dev.off()

############################################################################################
library(readxl)
x <- read_excel("~/Desktop/eNeuro/3results/ST2_GSEA_comparison.xlsx", 
                sheet = "GOCC")
names(x) <- gsub("CellAware", "corrected", names(x))
names(x) <- gsub("CellUnaware", "original", names(x))
x <- x[order(x$corrected_padj),]
x <- x[1:10, ]

x$path <- gsub("_", " ", x$path)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$corrected_padj)])
tiff("../4plots/GOCC_GSEA.tiff", units = "cm", width = 10, height = 12, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(corrected_padj)), data=x)+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(original_padj)), data=x)+
  geom_point( aes(x=path, y=-log10(original_padj), colour =original_NES), shape=16, size=4)+
  geom_point( aes(x=path, y=-log10(corrected_padj),  colour =corrected_NES), shape=18, size=4)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  ggtitle("GOCC")+ 
  theme(legend.position = "bottom",  legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed")+
  scale_color_gradient2(low = "#008b45e5", high = "#bb0021e5",  limit = c(-2.5,2.5)) + labs(colour = "NES") 

dev.off()