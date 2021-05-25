#############################################################
######### EWCE full network and central proteins  ###########
#############################################################

library("EWCE")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("GEOquery")

mydf <- read.delim2("../2resources/averageexpression.tsv", row.names = 1)
# mydf <- log2(mydf)
# mydf[mydf=="-Inf"] <- 0
l1 <- as.character(names(mydf))
l1 <- list(l1, l1)
gcd <- generate.celltype.data(mydf,
                              l1,
                              "SN")
load("../2resources/CellTypeData_SN.rda")

# CentralProteins 
cp<- read.table("../3results/centralproteins.txt", quote="\"", comment.char="")
cp <- cp[as.character(cp$V1) %in% rownames(ctd[[1]][["mean_exp"]]),]
cp <- as.character(cp)
# EWCE plot
ppilist <- list()
ppilist[["centralproteins"]] = EWCE::bootstrap.enrichment.test(sct_data = ctd, 
                                                               hits = as.character(cp), 
                                                               bg = rownames(ctd[[1]][["mean_exp"]]),
                                                               reps = 10000, 
                                                               annotLevel = 1, 
                                                               geneSizeControl = FALSE, 
                                                               genelistSpecies = "human", 
                                                               sctSpecies = "human")$results

ppilist[["centralproteins"]]$BHp <- p.adjust(ppilist[["centralproteins"]]$p, method = "BH")


## Full Network 
# CentralProteins 
# EWCE plot
ppiall <- read.csv("../3results/nodes_centrality.csv", row.names=1)
ppiall <- ppiall[as.character(ppiall$Row.names) %in% rownames(ctd[[1]][["mean_exp"]]),]
ppilist[["full"]] = EWCE::bootstrap.enrichment.test(sct_data = ctd, 
                                                    hits = as.character(ppiall$Row.names), 
                                                    bg = rownames(ctd[[1]][["mean_exp"]]),
                                                    reps = 10000, 
                                                    annotLevel = 1, 
                                                    geneSizeControl = FALSE, 
                                                    genelistSpecies = "human", 
                                                    sctSpecies = "human")$results

ppilist[["full"]]$BHp <- p.adjust(ppilist[["full"]]$p, method = "BH")

# Binding and plot
ppilist2 <- bind_rows(ppilist, .id="Cluster")
ppilist2$Notation <- ifelse(ppilist2$BHp<0.05,
                            paste0( signif(ppilist2$sd_from_mean,3), "\n (", signif(ppilist2$BHp,2), ")"),
                            "")
ppilist2$Cluster <- as.character(ppilist2$Cluster)

ppilist2$CellType <- as.character(ppilist2$CellType)
ppilist2$CellType <- ifelse(ppilist2$CellType=="Neuron", "Neurons", 
                            ifelse(ppilist2$CellType=="Endothelial", "Endothelial Cells",
                                   ifelse(ppilist2$CellType=="OPC", "OPCs", 
                                          ifelse(ppilist2$CellType=="ODC", "ODCs", 
                                                 ifelse(ppilist2$CellType=="Astrocyte", "Astrocytes", "Microglia")))))

ppilist2$BHp <- p.adjust(ppilist2$p, method = "BH")
ppilist2$Notation <- ifelse(ppilist2$BHp<0.05,
                            paste0( signif(ppilist2$sd_from_mean,3), "\n (", signif(ppilist2$BHp,2), ")"),
                            "")
ppilist2$Notation2 <- ifelse(ppilist2$BHp>0.05,
                             paste0( signif(ppilist2$sd_from_mean,3), "\n (", signif(ppilist2$BHp,2), ")"),
                             "")
ppilist2$Cluster = ifelse(ppilist2$Cluster == "full","Entire\nNetwork","Central\nProteins" )
ppilist2$Cluster <- factor(as.character(ppilist2$Cluster), 
                           levels = c("Entire\nNetwork",
                                      "Central\nProteins"))


library("ggplot2")
tiff("../4plots/PPI_EWCE.tiff", width = 6, height = 2.5, res = 300, units = "in")
ggplot(ppilist2, aes(CellType,Cluster)) +
  geom_tile(aes(fill=sd_from_mean), alpha=0.5) + 
  theme_classic()+
  geom_text(aes(label = Notation), lineheight = .7, size=3) +
  geom_text(aes(label = Notation2), lineheight = .7, size=3, col="darkgrey") +
  scale_fill_gradient2(low = "white", high = "red", mid = "white",
                       space = "Lab", 
                       name="SDs\nfrom\nmean")+
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank())+
  xlab("")+
  ylab("")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
dev.off()


library("openxlsx")
write.xlsx(ppilist2, file = "../3results/EDT3-2.xlsx")


#############################################################
############ GSEA beweenness-ordered PPI network  ###########
#############################################################

library("biomaRt")
library("fgsea")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
canonical <- gmtPathways('../2resources/canonical_all.gmt')
GOBP <- gmtPathways('../2resources/GOBP.gmt')
GOMF <- gmtPathways('../2resources/GOMF.gmt')
GOCC <- gmtPathways('../2resources/GOCC.gmt')
nodes_centrality <- read.csv("../5_PPI_API/cleaned/nodes_centrality.csv", row.names=1)
nodes_centrality <- nodes_centrality[order(nodes_centrality$betw, decreasing = T),]
nodes_centrality$ID <- 1:nrow(nodes_centrality)
nodes_centrality$betw <- as.numeric(as.character(nodes_centrality$betw))
nodes_centrality$Score<-(nodes_centrality$betw-min(nodes_centrality$betw))/(max(nodes_centrality$betw)-min(nodes_centrality$betw))
entrez <- getBM(attributes= c("hgnc_symbol", "entrezgene_id"),
                values = nodes_centrality$Row.names,
                mart = ensembl)
names(entrez) <- c("Gene", "ENTREZ" )
names(nodes_centrality)  <- c("Gene", "deg", "betw", "degrank", "betrank", "Score")
nodes_centrality <- merge(nodes_centrality, entrez, by="Gene")
nodes_centrality <- na.omit(nodes_centrality)
nodes_centrality <- nodes_centrality[order(-nodes_centrality$Score ), ]
ranks <- nodes_centrality$Score
names(ranks) <- nodes_centrality$ENTREZ

annot <- Hmisc::llist(GOCC, GOBP, GOMF, canonical)
gsea <- list()
for (i in names(annot)){
  gsea[[i]] <- fgsea(annot[[i]], ranks,  minSize = 25, maxSize=500, nperm=100000)
  gsea[[i]]$padj <- p.adjust(gsea[[i]]$pval, method = "BH")
}

library(openxlsx)
write.xlsx(gsea, file = paste0("../3results/GSEA_PPI.xlsx"))
###########################################################################
# PPI GSEA plots
library(readxl)
x <- read_excel("../3results/GSEA_PPI.xlsx", 
                sheet = "canonical")

library("ggplot2")
x <- x[order(x$padj),]
x <- x[1:5, ]
x$path <- gsub("_", " ", x$pathway)
x$path <- gsub("REACTOME", "REACTOME - ", x$path)
x$path <- gsub("KEGG", "KEGG - ", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$padj)])
tiff("../4plots/canonicalGSEA_PPI.tiff", units = "cm", width = 4.5, height = 3.5, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(padj)), data=x, size=0.2)+
  geom_point( aes(x=path, y=-log10(padj),  colour =NES), shape=18, size=2)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  theme(legend.position = "none",  
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 4),
        plot.title = element_text(size = 6),
        axis.line = element_line(colour = "black", size = 0.05),
        axis.ticks = element_line(colour = "black", size = 0.05))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.1)+
  labs(colour = "NES") + 
  scale_color_gradient2(low = "royalblue", high = "royalblue",  limit = c(-2.5,2.5))+
  ggtitle("Canonical data set")
dev.off()

############################################################################################
library(readxl)
x <- read_excel("../3results/GSEA_PPI.xlsx", 
                sheet = "GOCC")

library("ggplot2")
x <- x[order(x$padj),]
x <- x[1:5, ]
x$path <- gsub("_", " ", x$pathway)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$padj)])
tiff("../4plots/GOCCGSEA_PPI.tiff", units = "cm", width = 4.5, height = 3.5, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(padj)), data=x, size=0.2)+
  geom_point( aes(x=path, y=-log10(padj),  colour =NES), shape=18, size=2)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  theme(legend.position = "none",  
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 4),
        plot.title = element_text(size = 6),
        axis.line = element_line(colour = "black", size = 0.05),
        axis.ticks = element_line(colour = "black", size = 0.05))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.1)+
  labs(colour = "NES") + 
  scale_color_gradient2(low = "royalblue", high = "royalblue",  limit = c(-2.5,2.5))+
  ggtitle("GOCC")
dev.off()

############################################################################################
library(readxl)
x <- read_excel("../3results/GSEA_PPI.xlsx", 
                sheet = "GOMF")

library("ggplot2")
x <- x[order(x$padj),]
x <- x[1:5, ]
x$path <- gsub("_", " ", x$pathway)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$padj)])
tiff("../4plots/GOMFGSEA_PPI.tiff", units = "cm", width = 4.5, height = 3.5, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(padj)), data=x, size=0.2)+
  geom_point( aes(x=path, y=-log10(padj),  colour =NES), shape=18, size=2)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  theme(legend.position = "none",  
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 4),
        plot.title = element_text(size = 6),
        axis.line = element_line(colour = "black", size = 0.05),
        axis.ticks = element_line(colour = "black", size = 0.05))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.1)+
  labs(colour = "NES") + 
  scale_color_gradient2(low = "royalblue", high = "royalblue",  limit = c(-2.5,2.5))+
  ggtitle("GOMF")
dev.off()
############################################################################################
library(readxl)
x <- read_excel("../3results/GSEA_PPI.xlsx", 
                sheet = "GOBP")

library("ggplot2")
x <- x[order(x$padj),]
x <- x[1:5, ]
x$path <- gsub("_", " ", x$pathway)
x$path <- gsub("GO ", "", x$path)
x$path <- factor(x$path, levels = x$path[order(-x$padj)])
tiff("../4plots/GOBPGSEA_PPI.tiff", units = "cm", width = 4.5, height = 3.5, res = 300)
ggplot(x) +
  coord_flip()+
  geom_segment(aes(x = path, 
                   xend = path, 
                   y= 0,
                   yend = -log10(padj)), data=x, size=0.2)+
  geom_point( aes(x=path, y=-log10(padj),  colour =NES), shape=18, size=2)+
  ylab("-log10(BHp)") +
  xlab("")+
  theme_classic() +
  theme(legend.position = "none",  
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 4),
        plot.title = element_text(size = 6),
        axis.line = element_line(colour = "black", size = 0.05),
        axis.ticks = element_line(colour = "black", size = 0.05))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.1)+
  labs(colour = "NES") + 
  scale_color_gradient2(low = "royalblue", high = "royalblue",  limit = c(-2.5,2.5))+
  ggtitle("GOBP")
dev.off()


#############################################################
####### FET PPI/Central Proteins with PD/GWAS genes #########
#############################################################
nalls2019_gwas_catalog <- read.delim("../../2resources/nalls2019_gwas_catalog.tsv")
nalls2019_gwas_catalog <- unique(nalls2019_gwas_catalog$REPORTED.GENE.S.)
EDT4_1_xlsx <- read_excel("../3results/EDT4-1.xlsx")
EDT4_1_xlsx <- EDT4_1_xlsx$Node

library(GeneOverlap)
# Number of protein coding genes https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Homo_sapiens/109/
# GWAS proximal vs whole network 
go.obj <- newGeneOverlap(EDT4_1_xlsx,
                         nalls2019_gwas_catalog,
                         genome.size = 20203
)
go.obj <- testGeneOverlap(go.obj)
go.obj


nalls2019_gwas_catalog <- read.delim("../../2resources/nalls2019_gwas_catalog.tsv")
EDT4_1_xlsx <- read_excel("~/Desktop/eNeuro/manuscript/EDT4-1.xlsx.xlsx")
background <- EDT4_1_xlsx$Node
cp <- EDT4_1_xlsx$Central_Proteins

nalls2019_gwas_catalog <- nalls2019_gwas_catalog[nalls2019_gwas_catalog$REPORTED.GENE.S. %in% 
                                                   background,]$REPORTED.GENE.S.


# GWAS proximal vs Central proteins 
go.obj <- newGeneOverlap(cp,
                         intersect(nalls2019_gwas_catalog, EDT4_1_xlsx$Node),
                         genome.size = 5615)
go.obj <- testGeneOverlap(go.obj)
go.obj

########################################################################
PDg = c("LRRK2", "SNCA", "PARK7", "VPS35",  "PRKN", "PINK1")
# PD genes vs whole network  
go.obj <- newGeneOverlap(EDT4_1_xlsx$Node,
                         PDg,
                         genome.size=20203
)
go.obj <- testGeneOverlap(go.obj)
go.obj

# PD genes vs central proteins  
go.obj <- newGeneOverlap(cp,
                         intersect(EDT4_1_xlsx$Node, PDg),
                         genome.size=5615)
go.obj <- testGeneOverlap(go.obj)
go.obj


