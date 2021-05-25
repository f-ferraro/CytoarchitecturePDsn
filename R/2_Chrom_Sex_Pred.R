#####################################################################
######################### Gender imputation #########################
#####################################################################
library(tidyverse)
library(massiR)
library(biomaRt)
mart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
filters <- listFilters(mart)
attributes <- listAttributes(mart)
# Retrieve all genes
gene.attributes <-getBM(mart=mart, values=TRUE,
                        attributes= c("hgnc_symbol", "entrezgene_id",
                                      "chromosome_name", "start_position",
                                      "end_position", "strand"))
# Remove multiple mapping
unique.probe <- subset(gene.attributes, subset=!duplicated(gene.attributes[,1]))
# Keep only y 
y.unique <- subset(unique.probe, subset=unique.probe$chromosome_name == "Y")
# Make in the format for massiR
yprobes <- data.frame(row.names=y.unique$hgnc_symbol)
rm(list=setdiff(ls(), "yprobes"))

metadataall <- read.csv("../2resources/metadataall.csv", row.names = 1)
metadataall$GEO_series <- as.factor(metadataall$GEO_series)
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
  datExp[[paste(n, sep="")]] <- as.data.frame(mt)
  rm(mt, meta)
}

sample.results <- list()

# massiR functions
for (i in studies){
datExp1 <- datExp[[i]]
yc2 <- yprobes[rownames(yprobes) %in% rownames(datExp1),]
massi.y.out <- massi_y(datExp1, yc2)
massi_y_plot(massi.y.out)
massi.select.out <- massi_select(datExp1, yc2, threshold=4)
results <- massi_cluster(massi.select.out)
sample.results[[i]] <- data.frame(results[[2]])
}

rm(list=setdiff(ls(), "sample.results"))
metadataall <- read.csv("../2resources/metadataall.csv")

sample.results <- bind_rows(sample.results)

metadataall <- dplyr::select(metadataall, "ID", "GEO_series", "Gender")
names(metadataall) <- c("ID", "GEO_series", "Annot_Gender")
sample.results <- dplyr::select(sample.results, "ID", "sex")
names(sample.results) <- c("ID", "Predicted_Gender")

#
imputed <- merge(metadataall, sample.results, by="ID", all="TRUE")
imputed[imputed=="N/A"] <- NA
imputed$Predicted_Gender <- ifelse(imputed$Predicted_Gender=="female","F", "M")
imputed$Count <- ifelse(imputed$Annot_Gender == imputed$Predicted_Gender, 1, 0)
imputed2 <- na.omit(imputed)
sum(imputed2$Count, na.rm = T)/length(((imputed2$Count)))


metadataall <- read.csv("../2resources/metadataall.csv", row.names = 1)
metadataall[metadataall=="N/A"] <- NA

imputed$Annot_Gender <- imputed$Count <- NULL

metadataall <- merge(metadataall, imputed, by=c("ID", "GEO_series"), all=T)
metadataall$sumGender <- metadataall$Predicted_Gender
metadataall$sumGender <- ifelse(is.na(metadataall$sumGender),
                                       as.character(metadataall$Gender),
                                       as.character(metadataall$sumGender))

metadataall$X.1 <- metadataall$X <- NULL 
write.csv(metadataall, "../2resources/metadatall_predsex.csv", row.names = F)

# Counts
x <- metadataall %>% 
  group_by(GEO_series, sumGender) %>%
  summarise(count = n())

x <- pivot_wider(x, names_from = "sumGender", values_from = "count")



