#####################################################################
################### Marker selection from GSE140231 #################
#####################################################################

library("BRETIGEA")
library("ggsci")
library("gplots")
library("metafor")
library("tidyverse")
library("ggpubr")
library("ggplot2")
library("TOAST")
library("lmerTest")
library("broom")
library("broom.mixed")

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
SNC <- read.csv("../2resources/SNC.markers")
SNC$cluster <- gsub("GABA", "Neuron", SNC$cluster)
SNC$cluster <- gsub("DaNs", "Neuron", SNC$cluster)
marker <- SNC
marker <- dplyr::select(marker, "gene", "cluster", "avg_logFC")
marker <- marker[!duplicated(marker$gene), ]
names(marker) <- c("markers", "cell", "FC")

mrk <- marker
mrk <- mrk[order(-mrk$FC),] 

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

tiff("../4plots/SF1b_sc_VrProp.tiff", width = 10, height = 8, units = "in", res = 150)
ggplot(effect[as.numeric(as.character(effect$nmark))<150,], aes(as.numeric(as.character(nmark)),
                                                                as.numeric(as.character(pvar)),
                                                                color=Study))+
  geom_line()+
  facet_wrap(~Cell,ncol=3, nrow=2)+
  scale_x_continuous(breaks = seq(0,150, 10))+
  geom_hline(yintercept=0.35, linetype="dashed")+
  ggsci::scale_colour_npg()+
  # geom_vline(xintercept=10, linetype="dashed")+
  xlab("Number of markers included")+
  ylab("Proportion of total variance explained")+
  # theme(axis.text=element_text(size=5))+
  theme_classic(base_size = 10)
dev.off()
# save.image("../3results/variance_sc.RData")
load("../3results/variance_sc.RData")

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


#Top100
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
  mt <- exp(mt)
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
    filt2$pd <- rowSums(filt2 < 0) 
    for (it in 1:100){
      if (max(filt2$pd) >0) {
        filt2 <- filt2[filt2$pd < max(filt2$pd),]
        filt2 <- filt2[, names(filt2) %in% rownames(filt2)]
        filt2$pd <- rowSums(filt2 < 0)}}
    
    filt <- filt2[filt2$pd < 50, ] 
    keeping <- c(keeping, rownames(filt))
  }
  celltypes3[[study]]$mk <- rownames(celltypes3[[study]]) 
  celltypes3[[study]] <- subset(celltypes3[[study]], celltypes3[[study]]$mk %in% keeping)
  
}

for (n in as.character(studies)){
  meta <- subset(metadataall, metadataall$GEO_series == n)
  mt <- subset(totalmatrix, rownames(totalmatrix) %in% meta$ID)
  # mt <- scale(as.matrix(mt))
  mt <- as.data.frame(t(mt))
  mt <- na.omit(mt)
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.character))  
  mt[,1:ncol(mt)] <- as.matrix(sapply(na.omit(mt[,1:ncol(mt)]), as.numeric))  
  metastudy[[paste(n, sep="")]] <- meta
  mt <- exp(mt)
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  datExp[[paste(n, sep="")]] <- datExp[[paste(n, sep="")]][rownames(datExp[[paste(n, sep="")]]) %in% celltypes3[[n]]$mk,]
  rm(mt, meta)
}


datCor <- lapply(datExp, function(x)
  x <- cor(t(x),t(x)))


#############################################################################################
# Top20
ct <- celltypes3
ct <- lapply(ct, function(x){ x <- split(x, as.factor(x$cell));
x <- lapply(x, function(y){
  y <- as.character(rownames(y));
  y <- y[1:20];
  return(y)});
return(x)})

saveRDS(ct,"../3results/selectedMarkers_sc.rds")

ct2 <- lapply(ct, bind_cols)
ct2 <- lapply(ct2, function(x){
  x<- as.data.frame(x)
  x <- pivot_longer(x, -0, names_to="Cell", values_to="mkr")
  x<- as.data.frame(x)
  return(x)
})
#######################################
for (n in names(datCor)){
  
  datCor[[n]] <- as.data.frame(datCor[[n]])  
  datCor[[n]] <- datCor[[n]][rownames(datCor[[n]]) %in% ct2[[n]]$mkr,
                             names(datCor[[n]]) %in% ct2[[n]]$mkr]
  
}


library(pheatmap)
pdf("../4plots/CorMarkers.pdf", width = 7, height = 8)
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
