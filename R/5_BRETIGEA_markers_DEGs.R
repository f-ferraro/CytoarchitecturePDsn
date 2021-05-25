#####################################################################
################### Marker selection from BRETIGEA ##################
#####################################################################

library(BRETIGEA)
library(ggsci)
library(gplots)
library(metafor)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(TOAST)
library(lmerTest)
library(broom)
library(broom.mixed)

metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)

studies <- unique((as.character(metadataall$GEO_series)))

totalmatrix <- as.data.frame(t(totalmatrix))

metastudy <- list()
datExp <- list()

for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt <- exp(mt)
  metastudy[[paste(n, sep="")]] <- meta
  datExp[[paste(n, sep="")]] <- mt
  rm(mt, meta)
}

#
# ###############################################
# # Collapse to the major cell types
marker <- BRETIGEA::markers_df_human_brain

mrk <- marker
# Order by 1PC
for (type in unique(mrk$cell)){
  pcares <- list()
  effect <- list()
  
  ast <- mrk[mrk$cell == paste0(type), ]
  
  for (nmarker in 1:length(ast$markers)){
    
    effect[[nmarker]]<- list()
    subdatExp <- list()
    
    for (studies in names(datExp)){
      x <- subset(datExp[[studies]], rownames(datExp[[studies]]) %in% ast$markers)
      x <- x[is.finite(rowSums(x)),]
      x <- t(x)
      x<- x[ , which(apply(x, 2, var) != 0)]
      x <- t(x)
      pcares[[studies]] <- prcomp(t(x), scale. = TRUE)
      pcares[[studies]] <- as.data.frame(pcares[[studies]]$rotation)
      pcares[[studies]] <- pcares[[studies]][order(-abs(pcares[[studies]]$PC1)), ]
      pcares[[studies]] <- as.data.frame(pcares[[studies]])
      pcares[[studies]] <- pcares[[studies]][1:nmarker,]
      
      subdatExp[[studies]] <- subset(datExp[[studies]], rownames(datExp[[studies]]) %in% ast$markers)
      subdatExp[[studies]] <- subset(subdatExp[[studies]],
                                     rownames(subdatExp[[studies]]) %in% rownames(pcares[[studies]]))
      subdatExp[[studies]] <- t(subdatExp[[studies]])
      subdatExp[[studies]]<- subdatExp[[studies]][ , which(apply(subdatExp[[studies]], 2, var) != 0)]
      
      subdatExp[[studies]] <- prcomp((subdatExp[[studies]]), scale = TRUE)
    }
    
    for (studies in names(datExp)){
      effect[[nmarker]][[studies]]<- list()

      pca <- subdatExp[[studies]]
      PoV <- (pca$sdev^2/sum(pca$sdev^2))[1]
      
      effect[[nmarker]][[studies]] <- data.frame(pvar = PoV,
                                                 nmark = paste(nmarker))
    }
  }
  
  assign(paste0("pcares_", type), pcares)
  
  effect <- lapply(effect, function(x)
    x <- bind_rows(x, .id="Study"))
  effect <- bind_rows(effect)
  effect$Cell <- "cell"
  
  assign(paste0("effect_", type), effect)
}

# ## Plot to explore the proportion of variance per number of marker selected
Pattern1<-grep("effect_",names(.GlobalEnv),value=TRUE)
Pattern1_list<-do.call("list",mget(Pattern1))

effect <- bind_rows(Pattern1_list, .id="Cell")
effect$Cell <- gsub("effect_", "", effect$Cell)

tiff("../4plots/Br_VrProp.tiff", width = 13, height = 7, units = "in", res = 150)
ggplot(effect[as.numeric(as.character(effect$nmark))<150,], aes(as.numeric(as.character(nmark)),
                                                                as.numeric(as.character(pvar)),
                                                                color=Study))+
  geom_line()+
  facet_wrap(~Cell,ncol=3, nrow=2)+
  scale_x_continuous(breaks = seq(0,150, 10))+
  # geom_hline(yintercept=0.35, linetype="dashed")+
  ggsci::scale_colour_npg()+
  xlab("Number of markers included")+
  ylab("Proportion of total variance explained")+
  # theme(axis.text=element_text(size=5))+
  theme_classic(base_size = 10)
dev.off()
# save.image("../3results/varianceBR.RData")
load("../3results/varianceBR.RData")
## Marker quality control step
# To further ensure the marker quality, we produce the correlation matrix between the markers.
# We do not expect to see a perfect clustering between the cell types but still a good percentage of them should. 

celltypes<-grep("pcares_",names(.GlobalEnv),value=TRUE)
celltypes<-do.call("list",mget(celltypes))

rm(list=setdiff(ls(), "celltypes"))

celltypes <- lapply(celltypes, function(x) lapply(x, function(y)
  y<- rownames(y)))

celltypes <- lapply(celltypes, function(x) lapply(x, as.data.frame))
library(data.table)

celltypes <- lapply(celltypes, function(x) rbindlist(x, idcol = TRUE))
celltypes <- rbindlist(celltypes, idcol = TRUE)
names(celltypes) <- c("cell", "GEO_series", "markers")

celltypes <- split(celltypes, as.factor(celltypes$GEO_series))
celltypes <- lapply(celltypes, function(x) {x$GEO_series <- NULL; return(x)})


# Select top 100
celltypes1 <- lapply(celltypes, function(c){
  c <- split(c, as.factor(c$cell));
  c <- lapply(c, function(z){
    z <- as.data.frame(z);          
    z <- as.character(z$marker);
    z <- z[1:100] # Extraction value
    return(z)});
  return(c)})

celltypes1 <- lapply(celltypes1, function(x){
  names(x) <- gsub("pcares_", "", names(x));
  return(x)
})

cellypes2 <- lapply(celltypes1, function(s){
  s <- bind_cols(s)
  s <- pivot_longer(s, -0, names_to="cell", values_to="Value")
  s <- s$Value })

metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix <- as.data.frame(totalmatrix)
totalmatrix$Gene <- rownames(totalmatrix)

studies <- unique((as.character(metadataall$GEO_series)))

totalmatrix <- (t(totalmatrix))


metastudy <- list()
datExp <- list()


for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  datExp[[paste(n, sep="")]] <- datExp[[paste(n, sep="")]][rownames(datExp[[paste(n, sep="")]]) %in% cellypes2[[n]],]
  rm(mt, meta)
}


datCor <- lapply(datExp, function(x)
  x <- cor(t(x),t(x)))

celltypes3 <- lapply(celltypes1, function(s){
  s <- as.data.frame(s)
  s <- bind_cols(s)
  s <- pivot_longer(s, -0, names_to="cell", values_to="Value")
  s <- as.data.frame(s)
  s<- s[!duplicated(s$Value),]
  rownames(s) <- s$Value
  s$Value <- NULL;
  return(s)
})


# Filter for marker not in the same direction 
for (study in studies){
  cml <- as.data.frame(celltypes3[[study]]) 
  keeping <- as.character()
  for (cm in unique(cml$cell)){
    filt <- rownames(subset(cml, cml$cell %in% paste(cm)))
    df <- as.data.frame(datCor[[study]])
    filt2 <- df[ rownames(df) %in% filt, names(df) %in% filt ]
    filt2$pd <- rowSums(filt2 < 0) # We count the number of negatives in the correalation
    
    for (it in 1:100){
      if (max(filt2$pd) >0) {
        filt2 <- filt2[filt2$pd < max(filt2$pd),]
        filt2 <- filt2[, names(filt2) %in% rownames(filt2)]
        filt2$pd <- rowSums(filt2 < 0)}}
    
    filt <- filt2[filt2$pd < 50, ] # We filter out genes that show more than half negative correlations
    keeping <- c(keeping, rownames(filt))
  }
  celltypes3[[study]]$mk <- rownames(celltypes3[[study]]) 
  celltypes3[[study]] <- subset(celltypes3[[study]], celltypes3[[study]]$mk %in% keeping)
  
}

# Update correlation matrix 
for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  datExp[[paste(n, sep="")]] <- datExp[[paste(n, sep="")]][rownames(datExp[[paste(n, sep="")]]) %in% celltypes3[[n]]$mk,]
  rm(mt, meta)
}


#############################################################################################
# Top20  
ct <- celltypes3
ct <- lapply(ct, function(x){ x <- split(x, as.factor(x$cell));
x <- lapply(x, function(y){
  y <- as.character(rownames(y));
  y <- y[1:20];
  return(y)});
return(x)})

saveRDS(ct,"../3results/selectedMarkers.rds")
ct2 <- lapply(ct, bind_cols)
ct2 <- lapply(ct2, function(x){
  x<- as.data.frame(x)
  x <- pivot_longer(x, -0, names_to="Cell", values_to="mkr")
  x<- as.data.frame(x)
  return(x)
})
#######################################
# Update correlation matrix and plot
for (n in names(datCor)){
  
  datCor[[n]] <- as.data.frame(datCor[[n]])  
  datCor[[n]] <- datCor[[n]][rownames(datCor[[n]]) %in% ct2[[n]]$mkr,
                             names(datCor[[n]]) %in% ct2[[n]]$mkr]
  
}


library(pheatmap)
pdf("../4plots/SF2a_Br_CorMarkers.pdf", width = 7, height = 8)
for (study in names(datCor)){
  celltypes3[[study]]$mk <- NULL
  annotationb <- as.data.frame(celltypes3[[study]]) #
  # rownames(annotationb) <- rownames(annotation)
  names(annotationb) <- "Cell"
  
  pheatmap(datCor[[study]],
           annotation_col = annotationb,
           cluster_rows = T,
           cluster_cols = T,
           show_rownames = F,
           show_colnames = F, 
           main = study)
}
dev.off()
rm(list=ls())

#############################################################################################
################# Correlation among sn_markers and br_markers by study ######################
#############################################################################################
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix <- as.data.frame(totalmatrix)
totalmatrix$Gene <- rownames(totalmatrix)

studies <- unique((as.character(metadataall$GEO_series)))

totalmatrix <- (t(totalmatrix))


metastudy <- list()
datExp <- list()

# BR markers cor plots
ct <- readRDS("~/Scrivania/github/3results/selectedMarkers.rds")
ct2 <- lapply(ct, bind_cols)
ct2 <- lapply(ct2, function(x){
  x<- as.data.frame(x)
  x <- pivot_longer(x, -0, names_to="Cell", values_to="mkr")
  x<- as.data.frame(x)
  return(x)
})



for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  datExp[[paste(n, sep="")]] <- datExp[[paste(n, sep="")]][rownames(datExp[[paste(n, sep="")]]) %in% ct2[[n]]$mkr,]
  rm(mt, meta)
}


datCor <- lapply(datExp, function(x)
  x <- cor(t(x),t(x)))

bratigea <- list()

for (study in names(ct)){
  annotationb <- as.data.frame(ct2[[study]]) #
  rownames(annotationb) <- annotationb$mkr
  annotationb$mkr <- NULL
  names(annotationb) <- "Cell"
  
  tiff(paste0(study, "BRATIGEA.tiff"), units = "in", res = 300, height = 8, width = 8)   
  pheatmap(datCor[[study]],
           annotation_col = annotationb,
           cluster_rows = T,
           cluster_cols = T,
           show_rownames = F,
           show_colnames = F, 
           main = paste0("BRATIGEA\nmarkers ",study))
  dev.off()
}

rm(list=setdiff(ls(), "bratigea"))

############################################################################################################
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix <- as.data.frame(totalmatrix)
totalmatrix$Gene <- rownames(totalmatrix)

studies <- unique((as.character(metadataall$GEO_series)))

totalmatrix <- (t(totalmatrix))


metastudy <- list()
datExp <- list()

ct <- readRDS("../3results/selectedMarkers_sc.rds")

ct2 <- lapply(ct, bind_cols)
ct2 <- lapply(ct2, function(x){
  x<- as.data.frame(x)
  x <- pivot_longer(x, -0, names_to="Cell", values_to="mkr")
  x<- as.data.frame(x)
  return(x)
})



for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  datExp[[paste(n, sep="")]] <- datExp[[paste(n, sep="")]][rownames(datExp[[paste(n, sep="")]]) %in% ct2[[n]]$mkr,]
  rm(mt, meta)
}


datCor <- lapply(datExp, function(x)
  x <- cor(t(x),t(x)))

snsc <- list()

for (study in names(ct)){
  annotationb <- as.data.frame(ct2[[study]]) #
  annotationb <- na.omit(annotationb)
  rownames(annotationb) <- annotationb$mkr
  annotationb$mkr <- NULL
  names(annotationb) <- "Cell"
  
  tiff(paste0(study, "snsc.tiff"), units = "in", res = 300, height = 8, width = 8)  
  pheatmap(datCor[[study]],
           annotation_col = annotationb,
           cluster_rows = T,
           cluster_cols = T,
           show_rownames = F,
           show_colnames = F, 
           main = paste0("snRNAseq\nmarkers ", study))
  dev.off()
}

rm(list=ls())

#####################################################################
################ Cellular proportions deconvolution #################
#####################################################################
celltypes <- readRDS("../3results/selectedMarkers.rds")
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix$Gene <- rownames(totalmatrix)


totalmatrix <- (t(totalmatrix))

metastudy <- list()
datExp <- list()
for (n in unique(metadataall$GEO_series)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[n]] <- meta
  mt <- exp(mt)
  ct <- as.character(unlist(celltypes[[n]]))
  datExp[[n]] <- as.data.frame(mt)
  rm(mt, meta)
}

resCB <- list()
for (n in names(datExp)){
  resCB[[n]] <- MDeconv(datExp[[n]], celltypes[[n]],
                        epsilon = 1e-3, verbose = FALSE)$H
  resCB[[n]] <- as.data.frame(resCB[[n]])
}

#####################################################
# Cellular proportion correlation 
rm(list=setdiff(ls(), c("resCB", "celltypes")))

cor.res <- dplyr::bind_cols(resCB)
cor.res <- cor(t(cor.res), t(cor.res), use="pairwise.complete.obs", method = "spearman")
# Pvalues
cor.res <- as.data.frame(cor.res)
cor.res$Var1 <- rownames(cor.res)

cor.res <- pivot_longer(cor.res, -"Var1", names_to="Var2", values_to="r")

cor.res$n <- cor.res$pval <- cor.res$est <- 0

cor.res2 <- as.matrix(cor.res)
db <- dplyr::bind_cols(resCB)
db <- as.data.frame(t(db))

for (i in 1:nrow(cor.res)){
  test <- cor.res2[i,]
  from <- test[1]
  to <-test[2]
  ct <- cor.test(db[, from],
                 db[, to],
                 use = "pairwise.complete.obs", method="spearman")
  cor.res[i,]$est <- ct[["estimate"]][["rho"]]
  # cor.res[i,]$n <- ct[["parameter"]][["df"]]
  cor.res[i,]$pval <- (ct[["p.value"]])
}

cor.res$est <- signif(cor.res$est, 2)
cor.res$Sing <- ifelse(cor.res$pval < 0.05, "*","")

cor.res$label <- paste0(cor.res$est,  
                        cor.res$Sing)

cor_refor <- pivot_wider(cor.res[, c(1,2, 8)], names_from = 2, values_from = 3)
cor_refor <- as.data.frame(cor_refor)
rownames(cor_refor) <- cor_refor$Var1
cor_refor$Var1 <- NULL

cor.res$Var1 <- gsub("ast", "Astrocytes", cor.res$Var1)
cor.res$Var1 <- gsub("neu", "Neurons", cor.res$Var1)
cor.res$Var1 <- gsub("opc", "OPCs", cor.res$Var1)
cor.res$Var1 <- gsub("end", "Endothelial Cells", cor.res$Var1)
cor.res$Var1 <- gsub("oli", "ODCs", cor.res$Var1)
cor.res$Var1 <- gsub("mic", "Microglia", cor.res$Var1)

cor.res$Var2 <- gsub("ast", "Astrocytes", cor.res$Var2)
cor.res$Var2 <- gsub("neu", "Neurons", cor.res$Var2)
cor.res$Var2 <- gsub("opc", "OPCs", cor.res$Var2)
cor.res$Var2 <- gsub("end", "Endothelial Cells", cor.res$Var2)
cor.res$Var2 <- gsub("oli", "ODCs", cor.res$Var2)
cor.res$Var2 <- gsub("mic", "Microglia", cor.res$Var2)

cor.res$Notation =  paste(cor.res$label,"\n(",
                          signif(cor.res$pval, 1), ")", sep = "");

tiff("../4plots/Figure1B.tiff", width = 7, height = 6, res = 150, units = "in")
ggplot(cor.res, aes(Var2, Var1)) +
  geom_tile(aes(fill=r), alpha=0.5) + 
  theme_classic()+
  geom_text(aes(label = Notation)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman's\nCorrelation")+
  theme(axis.line=element_blank())+
  xlab("")+
  ylab("Cell Type")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
dev.off()

#####################################################
## Linear model for differences in cell proportions

res2 <- bind_cols(resCB)
res2$cell <- rownames(res2)
res2 <- as.data.frame(t(res2))
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
metadataall$Gender <- metadataall$Predicted_Gender
metadataall$X <- metadataall$X.1 <- NULL
res2$ID <- rownames(res2)
res2 <- merge(res2, metadataall, by="ID")


rc <- split(res2, as.factor(res2$GEO_series))

for (i in names(rc)){
  rc[[i]][, c(2)] <- as.numeric(as.character(rc[[i]][, c(2)]))
  rc[[i]][, c(3)] <- as.numeric(as.character(rc[[i]][, c(3)]))
  rc[[i]][, c(4)] <- as.numeric(as.character(rc[[i]][, c(4)]))
  rc[[i]][, c(5)] <- as.numeric(as.character(rc[[i]][, c(5)]))
  rc[[i]][, c(6)] <- as.numeric(as.character(rc[[i]][, c(6)]))
  rc[[i]][, c(7)] <- as.numeric(as.character(rc[[i]][, c(7)]))
  rc[[i]][, c(11)] <- as.numeric(as.character(rc[[i]][, c(11)]))}

for (AGE_GEN in c("GSE20164","GSE20292", "GSE20333", "GSE8397")){
  rc[[AGE_GEN]] <- dplyr::select(rc[[AGE_GEN]], "ast", "end", "mic", "neu", "oli", "opc",
                          "ID", "Status", "Age", "Gender")}
for (AGE_GEN_BRAAK in c("GSE43490", "GSE42966")){
  rc[[AGE_GEN_BRAAK]] <- dplyr::select(rc[[AGE_GEN_BRAAK]],  "ast", "end", "mic", "neu", "oli", "opc",
                                "ID", "Status" , "Age", "Gender", "Braak"
  )}
for (BRAAK_GEN in c("GSE49036")){
  rc[[BRAAK_GEN]] <- dplyr::select(rc[[BRAAK_GEN]],  "ast", "end", "mic", "neu", "oli", "opc",
                        "ID", "Status", "Gender", "Braak")}
for (GEN in c("GSE7621", "GSE20163")){
  rc[[GEN]] <- dplyr::select(rc[[GEN]], "ast", "end", "mic", "neu", "oli", "opc",
                      "ID", "Status", "Gender")}



rc[["GSE43490"]]$Braak <- gsub("-", ".", rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- gsub("BR", "", rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- ifelse(is.na(rc[["GSE43490"]]$Braak ), 1, rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- as.numeric(as.character(rc[["GSE43490"]]$Braak))

rc[["GSE42966"]]$Braak <- gsub("-", ".", rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- gsub("BR", "", rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- gsub("N/A", 1, rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- ifelse(is.na(rc[["GSE42966"]]$Braak ), 1, rc[["GSE42966"]]$Braak)

rc[["GSE49036"]]$Braak <- gsub("-", ".", rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub("BR", "", rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub("N/A", 1, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub(3.4, 3.5, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub(5.6, 5.5, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- ifelse(is.na(rc[["GSE49036"]]$Braak ), 1, rc[["GSE49036"]]$Braak)


lm <- list()
for (i in names(rc)){
  lm[[i]] <- list()
  rc[[i]]$ID <- NULL
  for(lin in c( "ast", "end", "mic", "neu", "oli", "opc")){
    
    Formula <- formula(paste(lin, " ~ ", 
                             paste(names(rc[[i]][, 7:ncol(rc[[i]])]), collapse=" + ")))
    lm[[i]][[lin]] <- tidy(lm(Formula, rc[[i]]))
    lm[[i]][[lin]]$psign <- ifelse(lm[[i]][[lin]]$p.value < 0.05, "*","")
  }}


p <- list()
scaleFUN <- function(x) sprintf("%.2f", x)

for (study in names(rc)){
  for(lin in c( "ast", "end", "mic", "neu", "oli", "opc")){
    x <- lm[[study]][[lin]]
    x <- subset(x, term=="StatusPD")
    dt <- rc[[study]]
    dt <- dt[, c("Status", paste(lin))]
    names(dt) <- c("Status", "Cell")
    p[[paste0(study, "_", lin)]] <-   ggplot(dt, 
                                             aes(x=Status, y=Cell, fill=Status)) + 
      geom_violin()+
      geom_jitter(height = 0, width = 0.04, size=0.1)+
      theme(strip.text.x = element_text(size = 8, colour = "black"))+
      theme_classic()+
      ggsci::scale_fill_lancet(alpha = 0.3)+
      xlab(NULL)+
      ylab(NULL) +
      annotate("text", x = 1.5, y = 0.75*max(dt$Cell), size= 8,  label = paste(x$psign)) +
      theme(legend.position = "none", axis.text=element_text(size=6))+
      stat_summary(fun=mean, geom="point", shape=23)+
      theme(  axis.title.y = element_text(size = 8))+
      scale_y_continuous(labels=scaleFUN)
    
  }
}


p$GSE20163_ast <- p$GSE20163_ast + ggtitle("Astrocytes")
p$GSE20163_end <- p$GSE20163_end + ggtitle("Endothelial Cells")
p$GSE20163_mic <- p$GSE20163_mic + ggtitle("Microglia")
p$GSE20163_neu <- p$GSE20163_neu+ ggtitle("Neurons")
p$GSE20163_oli <- p$GSE20163_oli + ggtitle("ODCs")
p$GSE20163_opc <- p$GSE20163_opc + ggtitle("OPCs") 
p$GSE20163_ast<- p$GSE20163_ast + ylab("GSE20163")
p$GSE20164_ast <- p$GSE20164_ast + ylab("GSE20164")
p$GSE20292_ast <- p$GSE20292_ast + ylab("GSE20292")
p$GSE20333_ast <- p$GSE20333_ast + ylab("GSE20333")
p$GSE42966_ast <- p$GSE42966_ast + ylab("GSE42966")
p$GSE43490_ast <- p$GSE43490_ast + ylab("GSE43490")
p$GSE49036_ast <- p$GSE49036_ast + ylab("GSE49036")
p$GSE7621_ast <- p$GSE7621_ast + ylab("GSE7621")
p$GSE8397_ast <- p$GSE8397_ast + ylab("GSE8397")

library(egg)
tiff("../4plots/Figure1A.tiff", width=10, height = 7, units="in", res=190)
do.call('ggarrange',c(p, ncol = 6))
dev.off()


lm <- lapply(lm, bind_rows, .id="cell")
lm <- bind_rows(lm, .id="Study")
lm <- lm[lm$term=="StatusPD",]
lm$cell <- ifelse(lm$cell=="neu", "Neurons", lm$cell)
lm$cell <- ifelse(lm$cell=="oli", "ODCs", lm$cell)
lm$cell <- ifelse(lm$cell=="opc", "OPCs", lm$cell)
lm$cell <- ifelse(lm$cell=="end", "Endothelial Cells", lm$cell)
lm$cell <- ifelse(lm$cell=="mic", "Microglia", lm$cell)
lm$cell <- ifelse(lm$cell=="ast", "Astrocytes", lm$cell)

library(openxlsx)

write.xlsx(lm, file = "../3results/EDT1-1.xlsx")


######################################################################################################
# half-violin for the cell proportions in the single studies
celltypes <- readRDS("../3results/selectedMarkers.rds")
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix$Gene <- rownames(totalmatrix)


totalmatrix <- (t(totalmatrix))

metastudy <- list()
datExp <- list()
for (n in unique(metadataall$GEO_series)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[n]] <- meta
  mt <- exp(mt)
  ct <- as.character(unlist(celltypes[[n]]))
  datExp[[n]] <- as.data.frame(mt)
  rm(mt, meta)
}

resCB <- list()
for (n in names(datExp)){
  resCB[[n]] <- MDeconv(datExp[[n]], celltypes[[n]],
                        epsilon = 1e-3, verbose = FALSE)$H
  resCB[[n]] <- as.data.frame(resCB[[n]])
}


res2 <- bind_cols(resCB)
res2$cell <- rownames(res2)
res2 <- as.data.frame(t(res2))
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
metadataall$Gender <- metadataall$Predicted_Gender
metadataall$X <- metadataall$X.1 <- NULL
res2$ID <- rownames(res2)
res2 <- merge(res2, metadataall, by="ID")


rc <- split(res2, as.factor(res2$GEO_series))

for (i in names(rc)){
  rc[[i]][, c(2)] <- as.numeric(as.character(rc[[i]][, c(2)]))
  rc[[i]][, c(3)] <- as.numeric(as.character(rc[[i]][, c(3)]))
  rc[[i]][, c(4)] <- as.numeric(as.character(rc[[i]][, c(4)]))
  rc[[i]][, c(5)] <- as.numeric(as.character(rc[[i]][, c(5)]))
  rc[[i]][, c(6)] <- as.numeric(as.character(rc[[i]][, c(6)]))
  rc[[i]][, c(7)] <- as.numeric(as.character(rc[[i]][, c(7)]))
  rc[[i]][, c(11)] <- as.numeric(as.character(rc[[i]][, c(11)]))}

for (AGE_GEN in c("GSE20164","GSE20292", "GSE20333", "GSE8397")){
  rc[[AGE_GEN]] <- dplyr::select(rc[[AGE_GEN]], "ast", "end", "mic", "neu", "oli", "opc",
                                 "ID", "Status", "Age", "Gender")}
for (AGE_GEN_BRAAK in c("GSE43490", "GSE42966")){
  rc[[AGE_GEN_BRAAK]] <- dplyr::select(rc[[AGE_GEN_BRAAK]],  "ast", "end", "mic", "neu", "oli", "opc",
                                       "ID", "Status" , "Age", "Gender", "Braak"
  )}
for (BRAAK_GEN in c("GSE49036")){
  rc[[BRAAK_GEN]] <- dplyr::select(rc[[BRAAK_GEN]],  "ast", "end", "mic", "neu", "oli", "opc",
                                   "ID", "Status", "Gender", "Braak")}
for (GEN in c("GSE7621", "GSE20163")){
  rc[[GEN]] <- dplyr::select(rc[[GEN]], "ast", "end", "mic", "neu", "oli", "opc",
                             "ID", "Status", "Gender")}



rc[["GSE43490"]]$Braak <- gsub("-", ".", rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- gsub("BR", "", rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- ifelse(is.na(rc[["GSE43490"]]$Braak ), 1, rc[["GSE43490"]]$Braak)
rc[["GSE43490"]]$Braak <- as.numeric(as.character(rc[["GSE43490"]]$Braak))

rc[["GSE42966"]]$Braak <- gsub("-", ".", rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- gsub("BR", "", rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- gsub("N/A", 1, rc[["GSE42966"]]$Braak)
rc[["GSE42966"]]$Braak <- ifelse(is.na(rc[["GSE42966"]]$Braak ), 1, rc[["GSE42966"]]$Braak)

rc[["GSE49036"]]$Braak <- gsub("-", ".", rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub("BR", "", rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub("N/A", 1, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub(3.4, 3.5, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- gsub(5.6, 5.5, rc[["GSE49036"]]$Braak)
rc[["GSE49036"]]$Braak <- ifelse(is.na(rc[["GSE49036"]]$Braak ), 1, rc[["GSE49036"]]$Braak)


lm <- list()
for (i in names(rc)){
  lm[[i]] <- list()
  rc[[i]]$ID <- NULL
  for(lin in c( "ast", "end", "mic", "neu", "oli", "opc")){
    
    Formula <- formula(paste(lin, " ~ ", 
                             paste(names(rc[[i]][, 7:ncol(rc[[i]])]), collapse=" + ")))
    lm[[i]][[lin]] <- tidy(lm(Formula, rc[[i]]))
    lm[[i]][[lin]]$psign <- ifelse(lm[[i]][[lin]]$p.value < 0.05, "*","")
  }}

rm(list=setdiff(ls(), c("lm", "rc")))

lm2 = lapply(lm, bind_rows, .id="Cell")
lm2 = bind_rows(lm2, .id="GEO_Study")
lm2 = lm2[lm2$term =="StatusPD", ]
rc3 <- lapply(rc, function(x) x[, c("ast", "end", "mic", "neu", "oli", "opc", "Status")])
rc3 = bind_rows(rc3, .id="GEO_Study")
rc3 = pivot_longer(rc3, -c("GEO_Study","Status"), names_to = "Cell", values_to = "Value")

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

rc3$Cell <- gsub("ast", "Astrocytes", rc3$Cell)
rc3$Cell <- gsub("neu", "Neurons", rc3$Cell)
rc3$Cell <- gsub("mic", "Microglia", rc3$Cell)
rc3$Cell <- gsub("oli", "ODCs", rc3$Cell)
rc3$Cell <- gsub("opc", "OPCs", rc3$Cell)
rc3$Cell <- gsub("end", "Endothelial\nCells", rc3$Cell)

tiff(paste0("../4plots/ViolinCEllProp.tiff"), units="in", res=300, width=10, height = 6)
ggplot(rc3, aes(y=Value, x=1, fill=Status))+
  geom_split_violin(alpha=0.4)+
  facet_grid(rows=vars(Cell),cols=vars(GEO_Study), scales = "free")+
  geom_jitter(rc3[rc3$Status=="CTRL",], mapping =aes(x=0.90, fill=Status, col=Status),alpha=0.7, size=0.8,
              position = position_jitter(w = 0.1, h = 0))+
  geom_jitter(rc3[rc3$Status=="PD",], mapping =aes(x=1.10, fill=Status, col=Status),alpha=0.7, size=0.8, 
              position = position_jitter(w = 0.1, h = 0))+
  stat_summary(fun=median, geom="point", shape=23, size=2)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 4))+
  theme_classic()+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  # scale_fill_manual(values=c("royalblue", "coral1"))+
  # scale_color_manual(values=c("royalblue", "coral1"))+
  xlab("Data set")+
  ylab("Cell Type")
dev.off()


###########################################################################
# Meta-analysis of the estimates
res2 <- bind_cols(resCB)
res2$cell <- rownames(res2)
res2 <- pivot_longer(res2, -"cell", names_to = "ID", values_to = "Estimate")
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
res2 <- merge(res2, metadataall, by="ID")

# Random-effects model 
rmdf <- dplyr::select(res2, "ID", "cell", "Estimate", "GEO_series", "Status")
rmdf <- res2 %>%
  group_by(GEO_series, cell, Status) %>%
  summarise(Average = mean(Estimate),
            SD = sd(Estimate), 
            n=n())

#Function to reoganize for CTRL and PD for metafor 
x_ct <- rmdf%>% 
  filter(Status =="CTRL")
names(x_ct) <- c("GEO_series","Cell","Status","CTRLmean","CTRLsd","CTRLn" )
x_ct$Status <- NULL
x_pd <- x_pd <-  rmdf%>% 
  filter(Status =="PD")
names(x_pd) <- c("GEO_series","Cell","Status","PDmean","PDsd","PDn" )
x_pd$Status <- NULL
merx <- merge(x_ct, x_pd, by=c("GEO_series", "Cell"))
rm(x_ct)
rm(x_pd)

merx <- escalc(measure="SMD", m2i=CTRLmean, sd2i=CTRLsd, n2i=CTRLn,
               m1i=PDmean, sd1i=PDsd, n1i=PDn, data=merx)

nz <- split(merx, as.factor(merx$Cell))

effects <- lapply(nz, rma.uni, method="REML") #To solve convergence issue
df <- data.frame(t(sapply(effects,c)))
df <- as.data.frame(df)

selection <- dplyr::select(df, 0, b, beta, k, se, zval, pval, ci.lb, 
                    ci.ub, tau2, se.tau2, QE, QEp, I2, H2)
selection <- as.data.frame(selection)
selection$pBH <- NA 
selection$pBH <- p.adjust(selection$pval, "BH")


selection$b <- as.numeric(as.character(selection$b))

rownames(selection) <- gsub("ast", "Astrocytes", rownames(selection))
rownames(selection) <- gsub("neu", "Neurons", rownames(selection))
rownames(selection) <- gsub("opc", "OPCs", rownames(selection))
rownames(selection) <- gsub("end", "Endothelial Cells", rownames(selection))
rownames(selection) <- gsub("oli", "ODCs", rownames(selection))
rownames(selection) <- gsub("mic", "Microglia", rownames(selection))

selection$lab <- rownames(selection)

selection <- selection[order(selection$b),]
selection$lab <- factor(selection$lab, levels = as.character(selection$lab))

selection$text1 <- selection$lab
selection$text2 <- selection$lab

selection$text_NEU <- ifelse(selection$lab=="Neurons", "Neurons", "")
selection$text_ALL <- ifelse(selection$lab!="Neurons", as.character(selection$lab), "")

selection$text2 <- paste0(" (BHp=", signif(selection$pBH, 1),")")
selection$text3 <- paste0(selection$lab, "\n", selection$text2)

selection$se <- as.numeric(as.character(selection$se))
selection$ci.ub <- as.numeric(as.character(selection$ci.ub))
selection$ci.lb <- as.numeric(as.character(selection$ci.lb))
selection$b <- as.numeric(as.character(selection$b))
selection$errorbar <- ifelse(selection$b>0,selection$b+selection$se, selection$b-selection$se)


tiff("../4plots/Figure1C.tiff", width = 5, height = 3.5, units = 'in', res = 300)
ggplot(selection, aes(lab, b, fill=lab)) +
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=b, ymax=errorbar), size=0.2)+
  ggsci::scale_fill_npg(alpha = 0.5)+
  theme_classic()+
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none" )+
  coord_flip()+
  xlab("")+
  ylab("Effect Size")+
  geom_text(data=selection,aes(x=lab, y=0, label=text_NEU), size=3.2, hjust=1, vjust=0)+
  # geom_text(data=selection,aes(x=ifelse(lab %in% c("Neurons", "ODCs", "Endothelial Cells"), lab, NA), 
  #                              y=ifelse(errorbar>0,errorbar+.03,errorbar-.03), label="*"),
  #                               size=9)+
  geom_text(data=selection,aes(x=lab, y=0, label=text_ALL), size=3.2, hjust=0, vjust=0)+
  geom_text(data=selection,aes(x=ifelse(lab=="Neurons", lab, NA ), y=0, label=paste0("(BHp=",signif(pBH,2), ")")),
            hjust=1, size=2, vjust=1.5)+
  geom_text(data=selection,aes(x=ifelse(lab!="Neurons", lab, NA ), y=0, label=paste0("(BHp=",signif(pBH,2), ")")), 
            hjust=0, size=2, vjust=1.5)
dev.off()

library(openxlsx)

write.xlsx(selection, file = "../3results/cellprop_metafor.xlsx")

rm(list=setdiff(ls(), ""))


##########################################################################
####################### cell-proportions unaware DEGs ####################
########################################################################## 
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
metadataall$Sex <- metadataall$sumGender
metadataall$Predicted_Gender <- metadataall$sumGender <- metadataall$Gender <- NULL
metadataall$ID <- as.character(metadataall$ID)
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix$Gene <- rownames(totalmatrix)
totalmatrix <- pivot_longer(totalmatrix, -"Gene", names_to = "ID", values_to = "Value")
totalmatrix <- na.omit(totalmatrix)
lng <- merge(totalmatrix, metadataall, by="ID")
rm(totalmatrix)
lost <- read.csv("lost.csv")
lng <- lng[lng$Gene %in% lost$x, ]
#####################################################################################
# DEGs identification original
mixed <- singular <- list()

for (i in unique(lng$Gene)){
  print(i)
  slng <- lng[lng$Gene ==i,]
  slng$Status <- ifelse(slng$Status=="PD", 1, 0)
  slng$Status <- as.factor(slng$Status)
  mixed[[i]] <- list()
  
  Formula <- formula(Value ~ 1+ Status+Sex+ (1|GEO_series))
  lmix <- lmer(Formula, 
               data=slng, REML = FALSE)
  
  if(isSingular(lmix)){singular[[i]] <- i}
  out <- as.data.frame(tidy(lmix))
  out$Gene <- i
  mixed[[i]] <- out
  
} 

save.image("../3results/Original.RData")
mixed2 <- mixed
mixed2 <- bind_rows(mixed2, .id="Gene")
singfilt <- bind_rows(singular, .id="Gene")
mixed2 <- mixed2[!(mixed2$Gene %in% as.data.frame(t(singfilt))$V1), ]

mixed2 <- split(mixed2, as.factor(mixed2$term))
mixed2 <- mixed2$Status1
mixed2$pBH <- p.adjust(mixed2$p.value, method = "BH")
write.csv2(mixed2, "../3results/originalDEGs.csv")
mixed2_f <- subset(mixed2, mixed2$pBH < 0.05)
original <- mixed2_f

rm(list=setdiff(ls(), "original"))

rm(list=ls())



##########################################################################
######################## cell-proportions-aware DEGs #####################
########################################################################## 
# Deconvolve again
celltypes <- readRDS("../3results/selectedMarkers.rds")
metadataall <- read.csv("../2resources/metadatall_predsex.csv")
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix$Gene <- rownames(totalmatrix)

studies <- unique((as.character(metadataall$GEO_series)))

totalmatrix <- (t(totalmatrix))

metastudy <- list()
datExp <- list()

for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  mt <- exp(mt)
  ct <- as.character(unlist(celltypes[[n]]))
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  rm(mt, meta)
}

resCB <- list()
for (n in names(datExp)){
  resCB[[n]] <- MDeconv(datExp[[n]], celltypes[[n]],
                        epsilon = 1e-3, verbose = FALSE)$H
  resCB[[n]] <- as.data.frame(resCB[[n]])
}


rm(list=setdiff(ls(), c("resCB")))
resCB2 <- as.data.frame(t(bind_cols(resCB)))
resCB2$ID <- rownames(resCB2)
#######################################################################################
# Merging 
metadataall <- read.csv("../2resources/metadatall_predsex.csv")

metadataall$Sex <- metadataall$sumGender
metadataall$Predicted_Gender <- metadataall$sumGender <- metadataall$Gender <- NULL
metadataall$ID <- as.character(metadataall$ID)
totalmatrix <- read.csv("../2resources/totalmatrix.csv", row.names=1)
totalmatrix$Gene <- rownames(totalmatrix)
totalmatrix <- pivot_longer(totalmatrix, -"Gene", names_to = "ID", values_to = "Value")
totalmatrix <- na.omit(totalmatrix)
lng <- merge(totalmatrix, metadataall, by="ID")

lng <- merge(lng, resCB2, by="ID")
rm(totalmatrix, metadataall, resCB2)


lng$GEO_series <- as.character(lng$GEO_series)
# Remove genes from a single study
lng2 <- lng %>% 
  group_by(Gene, GEO_series) %>% 
  summarise(n=n()) 

lng3 <- lng2 %>% 
  group_by(Gene) %>% 
  summarise(n=n()) %>%
  filter(n>1)

resCB2 <- bind_cols(resCB)
resCB2 <- as.data.frame(t(resCB2))
resCB2$ID <- rownames(resCB2)
lng <- subset(lng, lng$Gene %in% lng2$Gene )


mixed <- singular <- list()
# DEGs identification original
for (i in unique(lng$Gene)){
  print(i)
  slng <- lng[lng$Gene ==i,]
  slng$Status <- ifelse(slng$Status=="PD", 1, 0)
  slng$Status <- as.factor(slng$Status)
  mixed[[i]] <- list()
  
  Formula <- formula(Value ~ 1+ Status+neu+oli+opc+Sex+ (1|GEO_series))
  lmix <- lmer(Formula, 
               data=slng, REML = FALSE)
  
  if(isSingular(lmix)){singular[[i]] <- i}
  out <- as.data.frame(tidy(lmix))
  out$Gene <- i
  mixed[[i]] <- out
  
} 
save.image("../3results/Corrected.RData")
mixed2 <- mixed
mixed2 <- bind_rows(mixed2, .id="Gene")
singfilt <- bind_rows(singular, .id="Gene")
mixed2 <- mixed2[!(mixed2$Gene %in% as.data.frame(t(singfilt))$V1), ]

mixed2 <- split(mixed2, as.factor(mixed2$term))
mixed2 <- mixed2$Status1
mixed2$pBH <- p.adjust(mixed2$p.value, method = "BH")
mixed2_f <- subset(mixed2, mixed2$pBH < 0.05)
write.csv2(mixed2, "../3results/correctedDEGs.csv")
correctedDEGs  <- subset(correctedDEGs , correctedDEGs $pBH < 0.05)
correctedDEGs  <- subset(correctedDEGs, abs(correctedDEGs$estimate) > log(1.2) )
write.csv2(correctedDEGs, "../3results/correctedDEGs.csv")
rm(list=setdiff(ls(), "resCB"))


###########################################################################
# Expression plots comparison before and after cell-proportion adjustments 
load("../3results/Corrected.RData")
  corrected <- lng
  rm(list=setdiff(ls(), c("corrected")))
  
correctedDEGs <- read.csv2("../3results/correctedDEGs.csv", row.names=1)
correctedDEGs <- correctedDEGs[abs(correctedDEGs$estimate)>log(1.2) & correctedDEGs$pBH<0.05, ]

corrected$Corrected <- ""
corrected$Original <- ""
  
corrected$Gene <- as.character(corrected$Gene)
corrected <- corrected[as.character(corrected$Gene) %in% correctedDEGs$Gene,]

for (i in unique(correctedDEGs$Gene)){
    print(i)
    corrected2 <- corrected[as.character(corrected$Gene) ==i,]
    
    corrected2$Status <- ifelse(corrected2$Status=="PD", 1, 0)
    corrected2$Sex <- ifelse(corrected2$Sex=="M", 1, 0)
    corrected2$Status <- as.factor(corrected2$Status)
    corrected2$Sex <- as.factor(corrected2$Sex)
    
    corrected2 <- split(corrected2, as.factor(corrected2$GEO_series))
    
    
    for (study in names(corrected2)){
    Formula <- formula(Value ~ neu+oli+opc+as.factor(Sex))
    lmix <- lm(Formula, 
                 data=corrected2[[study]])
    corrected$Corrected <-ifelse(corrected$Gene==i & corrected$GEO_series==study, resid(lmix), corrected$Corrected)
    
    Formula <- formula(Value ~ as.factor(Sex))
    lmix <- lm(Formula, 
               data=corrected2[[study]])
    corrected$Original <-ifelse(corrected$Gene==i & corrected$GEO_series==study, resid(lmix) ,corrected$Original)
    
    
    }
    
    } 
  
  
rm(list=setdiff(ls(), c("corrected", "c_sp", "o_sp")))
corrected$Original <- as.numeric(as.character(corrected$Original)) 
corrected$Corrected <- as.numeric(as.character(corrected$Corrected)) 
  
  
corrected <- corrected[!is.na(corrected$Corrected),]

library("tidyverse")
corrected <- corrected[, c("ID", "Gene", "GEO_series", "Status", "Original", "Corrected")]
lc <- pivot_longer(corrected, -c("ID", "Gene", "Status", "GEO_series"), values_to="Value", names_to = "Model")
lc <- na.omit(lc)
library("ggpubr")
  
# Function for split violins
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                             draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                               
                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                 stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                           1))
                                 quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                 aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                 aesthetics$alpha <- rep(1, nrow(quantiles))
                                 both <- cbind(quantiles, aesthetics)
                                 quantile_grob <- GeomPath$draw_panel(both, ...)
                                 ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                 ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                             })
  
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }
  
  
##################################################################################################
  
lc$Model <- as.character(lc$Model)
lc$Model <- factor(lc$Model, levels = c("Original", "Corrected"))
  
# lc$Value <- exp(lc$Value)
  
for(i in unique(corrected$Gene)){
  tiff(paste0("../4plots/singleexp/",i,".tiff"), units="in", res=300, 9, 3)
  plot(  
  ggplot(lc[as.character(lc$Gene) ==i,], aes(y=Value, x=Status, fill=Model))+ geom_split_violin(alpha=0.5)+
    facet_grid(cols=vars(GEO_series), scales = "free")+
      stat_summary(fun=mean, geom="point", shape=23, size=2)+
      ggtitle(paste(i))+
      theme(legend.position = "none",
            strip.text.x = element_text(size = 6),
            axis.text.x = element_text(size = 4),
            axis.text.y = element_text(size = 4))+
      theme_classic()+
     scale_fill_manual(values=c("coral1", "royalblue"))
  )
    dev.off()
}
