#############################################################
#################### Sample download ########################
#############################################################

library("GEOquery")
library("tidyverse")

dn <- getGEO(GEO='GSE42966', GSEMatrix = TRUE)
GSE42966 <- dn[["GSE42966_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE42966 <- cbind(dn[["GSE42966_series_matrix.txt.gz"]]@featureData@data$GENE, GSE42966)
GSE42966 <- as.data.frame(GSE42966)
rm(dn)
GSE42966 <- apply(GSE42966, 2, as.character)
GSE42966 <- apply(GSE42966, 2, as.numeric)
GSE42966 <- as.data.frame(GSE42966)
GSE42966 <- subset(GSE42966, !is.na(GSE42966$V1))

dn <- getGEO(GEO="GSE43490", GSEMatrix = TRUE)
GSE43490 <- dn[["GSE43490_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE43490 <- cbind(dn[["GSE43490_series_matrix.txt.gz"]]@featureData@data$GENE, GSE43490)
GSE43490 <- as.data.frame(GSE43490)
rm(dn)
GSE43490$GSM1294104 <- GSE43490$GSM1294105 <- GSE43490$GSM1294106 <-
  GSE43490$GSM1294107 <- GSE43490$GSM1294108 <- GSE43490$GSM1294109 <- 
  GSE43490$GSM1294110 <- GSE43490$GSM1294111 <- GSE43490$GSM1294112 <- 
  GSE43490$GSM1294113 <- GSE43490$GSM1294114 <- GSE43490$GSM1294115 <- 
  GSE43490$GSM1294116 <- GSE43490$GSM1294117 <- GSE43490$GSM1294131 <- 
  GSE43490$GSM1294132 <- GSE43490$GSM1294133 <- GSE43490$GSM1294134 <- 
  GSE43490$GSM1294135 <- GSE43490$GSM1294136 <- GSE43490$GSM1294137 <- 
  GSE43490$GSM1294138 <- GSE43490$GSM1294139 <- GSE43490$GSM1294140 <- 
  GSE43490$GSM1294141 <- GSE43490$GSM1294142 <- GSE43490$GSM1294143 <- 
  GSE43490$GSM1294144  <- NULL
GSE43490 <- apply(GSE43490, 2, as.character)
GSE43490 <- apply(GSE43490, 2, as.numeric)
GSE43490 <- as.data.frame(GSE43490)
GSE43490 <- subset(GSE43490, !is.na(GSE43490$V1))



dn <- getGEO(GEO="GSE20333", GSEMatrix = TRUE)
GSE20333 <- dn[["GSE20333_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE20333 <- cbind(dn[["GSE20333_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID	, GSE20333)
GSE20333 <- as.data.frame(GSE20333)
rm(dn)
GSE20333 <- apply(GSE20333, 2, as.character)
GSE20333 <- apply(GSE20333, 2, as.numeric)
GSE20333 <- as.data.frame(GSE20333)
GSE20333 <- subset(GSE20333, !is.na(GSE20333$V1))

dn <- getGEO(GEO="GSE8397", GSEMatrix = TRUE, getGPL = TRUE)
GSE8397 <- dn[["GSE8397-GPL96_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE8397 <- cbind(dn[["GSE8397-GPL96_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID, GSE8397)
rm(dn)
GSE8397<- as.data.frame(GSE8397)
GSE8397 <- dplyr::select(GSE8397, "V1","GSM208633", "GSM208634", "GSM208630", "GSM208645", 
                         "GSM208631", "GSM208632", "GSM208635", "GSM208638", "GSM208644", 'GSM208639',
                         "GSM208641", "GSM208643", "GSM208637", "GSM208642", "GSM208640", "GSM208636")
GSE8397 <- apply(GSE8397, 2, as.character)
GSE8397 <- apply(GSE8397, 2, as.numeric)
GSE8397 <- as.data.frame(GSE8397)
GSE8397 <- subset(GSE8397, !is.na(GSE8397$V1))


dn <- getGEO(GEO="GSE20292", GSEMatrix = TRUE)
GSE20292 <- dn[["GSE20292_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE20292 <- cbind(dn[["GSE20292_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID , GSE20292)
GSE20292 <- as.data.frame(GSE20292)
rm(dn)
GSE20292$GSM606624 <-GSE20292$GSM606625 <- GSE20292$GSM606626 <- NULL
GSE20292 <- apply(GSE20292, 2, as.character)
GSE20292 <- apply(GSE20292, 2, as.numeric)
GSE20292 <- as.data.frame(GSE20292)
GSE20292 <- subset(GSE20292, !is.na(GSE20292$V1))

dn <- getGEO(GEO="GSE20163", GSEMatrix = TRUE)
GSE20163 <- dn[['GSE20163_series_matrix.txt.gz']]@assayData[["exprs"]]
GSE20163 <- cbind(dn[["GSE20163_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID , GSE20163)
GSE20163 <-as.data.frame(GSE20163)
rm(dn)
GSE20163 <- apply(GSE20163, 2, as.character)
GSE20163 <- apply(GSE20163, 2, as.numeric)
GSE20163 <- as.data.frame(GSE20163)
GSE20163 <- subset(GSE20163, !is.na(GSE20163$V1))

dn <- getGEO(GEO="GSE20164", GSEMatrix = TRUE)
GSE20164 <- dn[["GSE20164_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE20164 <- cbind(dn[["GSE20164_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID , GSE20164)
GSE20164 <- as.data.frame(GSE20164)
rm(dn)
GSE20164 <- apply(GSE20164, 2, as.character)
GSE20164 <- apply(GSE20164, 2, as.numeric)
GSE20164 <- as.data.frame(GSE20164)
GSE20164 <- subset(GSE20164, !is.na(GSE20164$V1))

dn <- getGEO(GEO="GSE7621", GSEMatrix = TRUE)
GSE7621 <- dn[["GSE7621_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE7621 <- cbind(dn[["GSE7621_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID , GSE7621)
GSE7621 <- as.data.frame(GSE7621)
rm(dn)
GSE7621 <- apply(GSE7621, 2, as.character)
GSE7621 <- apply(GSE7621, 2, as.numeric)
GSE7621 <- as.data.frame(GSE7621)
GSE7621 <- subset(GSE7621, !is.na(GSE7621$V1))

dn <- getGEO(GEO="GSE49036", GSEMatrix = TRUE)
GSE49036 <- dn[["GSE49036_series_matrix.txt.gz"]]@assayData[["exprs"]]
GSE49036 <- cbind(dn[["GSE49036_series_matrix.txt.gz"]]@featureData@data$ENTREZ_GENE_ID , GSE49036)
rm(dn)
GSE49036 <- as.data.frame(GSE49036)
GSE49036$GSM1192699 <- NULL
GSE49036$GSM1192700 <- NULL
GSE49036$GSM1192701 <- NULL
GSE49036$GSM1192702 <- NULL
GSE49036$GSM1192703 <- NULL
GSE49036 <- apply(GSE49036, 2, as.character)
GSE49036 <- apply(GSE49036, 2, as.numeric)
GSE49036 <- as.data.frame(GSE49036)
GSE49036 <- subset(GSE49036, !is.na(GSE49036$V1))

#############################################################
################ Remove duplicated probes ###################
#############################################################
library("WGCNA")

rowGroup <- as.character(GSE20163$V1)
GSE20163$V1 <- NULL
rowID <- as.character(rownames(GSE20163))
GSE20163 <- sapply( GSE20163, as.numeric )
GSE20163 <- collapseRows(GSE20163, rowGroup, rowID, 
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE20163 <- as.data.frame(GSE20163[["datETcollapsed"]])


rowGroup <- as.character(GSE20164$V1)
GSE20164$V1 <- NULL
rowID <- as.character(rownames(GSE20164))
GSE20164 <- sapply( GSE20164, as.numeric )
GSE20164 <- collapseRows(GSE20164,rowGroup, rowID, 
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE20164 <- as.data.frame(GSE20164[["datETcollapsed"]])


rowGroup <- as.character(GSE20292$V1)
GSE20292$V1 <- NULL
rowID <- as.character(rownames(GSE20292))
GSE20292 <- sapply( GSE20292, as.numeric )
GSE20292 <- collapseRows(GSE20292, rowGroup, rowID, 
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE20292 <- as.data.frame(GSE20292[["datETcollapsed"]])


rowGroup <- as.character(GSE20333$V1)
GSE20333$V1 <- NULL
rowID <- as.character(rownames(GSE20333))
GSE20333 <- sapply( GSE20333, as.numeric )
GSE20333 <- WGCNA::collapseRows(GSE20333, rowGroup, rowID, 
                                method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE20333 <- as.data.frame(GSE20333[["datETcollapsed"]])

rowGroup <- as.character(GSE42966$V1)
rowID <- as.character(rownames(GSE42966))
GSE42966$V1 <- NULL
GSE42966 <- sapply( GSE42966, as.numeric )
GSE42966 <- collapseRows(GSE42966, rowGroup , rowID,
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE42966 <- as.data.frame(GSE42966[["datETcollapsed"]])

rowGroup <- as.character(GSE43490$V1)
GSE43490$V1 <- NULL
rowID <- as.character(rownames(GSE43490))
GSE43490 <- sapply( GSE43490, as.numeric )
GSE43490 <- collapseRows(GSE43490,  rowGroup, rowID,
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE43490 <- as.data.frame(GSE43490[["datETcollapsed"]])


rowGroup <- as.character(GSE49036$V1)
GSE49036$V1 <- NULL
rowID <- as.character(rownames(GSE49036))
GSE49036 <- sapply( GSE49036, as.numeric )
GSE49036 <- collapseRows(GSE49036,  rowGroup, rowID, 
                         method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE49036 <- as.data.frame(GSE49036[["datETcollapsed"]])


rowGroup <- as.character(GSE8397$V1)
GSE8397$V1 <- NULL
rowID <- as.character(rownames(GSE8397))
GSE8397 <- sapply( GSE8397, as.numeric )
GSE8397 <- collapseRows(GSE8397,  rowGroup, rowID,
                        method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE8397 <- as.data.frame(GSE8397[["datETcollapsed"]])


rowGroup <- as.character(GSE7621$V1)
GSE7621$V1 <- NULL
rowID <- as.character(rownames(GSE7621))
GSE7621 <- sapply( GSE7621, as.numeric )
GSE7621 <- collapseRows(GSE7621,  rowGroup, rowID,
                        method="maxRowVariance", connectivityBasedCollapsing=TRUE)
GSE7621 <- as.data.frame(GSE7621[["datETcollapsed"]])

#####################################################################
######################## log normalization ##########################
#####################################################################
GSE20163 <- log(GSE20163)
GSE20164 <- log(GSE20164)
GSE20292 <- log(GSE20292)
GSE20333 <- log(GSE20333)
GSE42966 <- log(GSE42966)
GSE43490 <- log(GSE43490)
GSE49036 <- log(GSE49036)
GSE7621 <- log(GSE7621)
GSE8397 <- log(GSE8397)

#############################################################
# Create the final matrix 

multimerge <- function (mylist) {
  unames <- unique(unlist(lapply(mylist, rownames)))
  n <- length(unames)
  out <- lapply(mylist, function(df) {
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    return(tmp)
  })
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(bigout)
}

mylist <- list(GSE20163, GSE20164, GSE20292,
               GSE42966, GSE43490, GSE49036, 
               GSE20333, GSE7621, GSE8397)

Mydata <- multimerge(mylist)

mydata <- as.data.frame(Mydata)
names(mydata) <- gsub("_", "", names(mydata)) 


rm(list=setdiff(ls(), "mydata"))   

#####################################################################
###################### Quantile normalization #######################
#####################################################################
mydata$V1 <- NULL
mydata <-limma::normalizeBetweenArrays(mydata, method = "quantile")
mydata <- as.data.frame(mydata)


rm(list=setdiff(ls(), "mydata"))   

#####################################################################
######################### ENSEMBL to Symbol #########################
#####################################################################
mydata$V1 <- rownames(mydata)

library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
entrez <- getBM(attributes= c("hgnc_symbol", "entrezgene_id"),
                values = mydata$V1,
                mart = ensembl)
names(entrez) <- c("Gene", "V1" )

mydata <- merge(mydata, entrez, by.x="V1")
mydata$V1  <- NULL
mydata <- subset(mydata, !is.na(Gene), )
rowGroup <- as.character(mydata$Gene)
rowID <- as.character(rownames(mydata))
mydata$Gene<- NULL
mydata <- sapply(mydata[, 1:ncol(mydata)], as.numeric)
mydata <- collapseRows(mydata,  rowGroup, rowID, 
                       method="maxRowVariance", connectivityBasedCollapsing=TRUE)
mydata <- as.data.frame(mydata[["datETcollapsed"]])

mydata <- as.data.frame(mydata)

write.csv(mydata, "../2resources/totalmatrix.csv")
