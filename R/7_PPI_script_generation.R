#############################################################
####### Harmonization of data from the data bases ###########
#############################################################

# With this script we want to construct the ppi network of the DEGs from the corrected lmm
library(biomaRt)
huri <- read.table("../5_PPI_API/DN/HuRi/HI-union.tsv")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensmgene <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id"),
                 mart = ensembl)

# Convert aliases
library("AnnotationDbi")
library("DBI")
library("org.Hs.eg.db")
dbCon <- org.Hs.eg_dbconn()
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
names(aliasSymbol) = gsub("_", "",names(aliasSymbol))

aliasSymbol <- as.data.frame(aliasSymbol)
aliasSymbol = aliasSymbol[,c(2,5)]
names(aliasSymbol) <- c("V1", "symbol" )
names(ensmgene) = c("V1", "ensG")

x = merge(aliasSymbol, ensmgene, by="V1")
x$V1 = NULL
x = x[,c(2,1)]

names(ensmgene) = c("symbol", "V1")
x = merge(huri, ensmgene, by="V1")
names(ensmgene) = c("symbol", "V2")
x = merge(x, ensmgene, by="V2")
rm(list=setdiff(ls(), "x"))
peptides <- read.table("../5_PPI_API/DN/StringAPI/peptides.txt", quote="\"", comment.char="")

names(x) = c("ens1", "ens2", "Prot1", "Prot2")
x = x[, c("Prot1", "Prot2")]


x2 <- x[(x$Prot1 %in% peptides$V1) | (x$Prot2 %in% peptides$V1), ]
list <- unique(c(as.character(x2$Prot1), as.character(x2$Prot2)))
huri <- x[(x$Prot1 %in% list) & (x$Prot2 %in% list), ]
huri$Source <- "HuRI"
names(huri) <- c("Prot1", "Prot2", "Source")
write.csv(huri, "../5_PPI_API/cleaned/huri.csv")

STRING <- read.delim("../5_PPI_API/DN/StringAPI/STRING.tsv")
STRING <- subset(STRING, STRING$ncbiTaxonId == "9606")
STRING <- dplyr::select(STRING, preferredName_A, preferredName_B, escore)
STRING$escore <- as.numeric(as.character(STRING$escore))
STRING <- subset(STRING, STRING$escore >0)
STRING$escore <- NULL
names(STRING) <- c("Prot1", "Prot2")
STRING$Source <- "STRING"
write.csv(STRING, "../5_PPI_API/cleaned/string.csv")

IID_PPIs <- read.delim("../5_PPI_API/DN/IID/PPIs.txt")
IID_PPIs$UniProt1 <- IID_PPIs$UniProt2 <- NULL
IID_PPIs <- subset(IID_PPIs, IID_PPIs$evidence.type == "exp")
IID_PPIs$evidence.type <- NULL
names(IID_PPIs) <- c("Prot1", "Prot2")
IID_PPIs$Source <- "IID_PPIs"
write.csv(IID_PPIs, "../5_PPI_API/cleaned/IID_PPIs.csv")
#
#
biogrid <- read.delim("../5_PPI_API/DN/BioGrid/biogrid2.csv", header=FALSE)
biogrid  <- dplyr::select(biogrid, V3,V4, V10)
biogrid  <- subset(biogrid, biogrid$V10 =="9606")
biogrid$V10 <- NULL
peptides <- read.table("../5_PPI_API/DN/StringAPI/peptides.txt", quote="\"", comment.char="")
names(biogrid) = c("Prot1", "Prot2")
x2 <- biogrid[(biogrid$Prot1 %in% peptides$V1) | (biogrid$Prot2 %in% peptides$V1), ]
list <- unique(c(as.character(x2$Prot1), as.character(x2$Prot2)))
biogrid <- biogrid[(biogrid$Prot1 %in% list) & (biogrid$Prot2 %in% list), ]
names(biogrid) <- c("Prot1", "Prot2")
biogrid$Source <- "biogrid"
write.csv(biogrid, "../5_PPI_API/cleaned/biogrid.csv")
#
#
intact <- read.delim("../5_PPI_API/DN/IntAct/intact.csv", header=FALSE, comment.char="#")
intact <- subset(intact, intact$V12 %in% c("psi-mi:MI:0914(association)",
                                           "psi-mi:MI:0407(direct interaction",
                                           "psi-mi:MI:0915(physical association)"))
intact <- subset(intact, intact$V10 == "taxid:9606(human)|taxid:9606(Homo sapiens)")
intact <- subset(intact, intact$V11 == "taxid:9606(human)|taxid:9606(Homo sapiens)")
intact <- dplyr::select(intact, V1, V2)
library("stringr")
intact$V1 <- word(intact$V1,2,sep = "\\:")
intact$V2 <- word(intact$V2,2,sep = "\\:")
intact$V1 <- word(intact$V1,1,sep = "\\-")
intact$V2 <- word(intact$V2,1,sep = "\\-")
library("biomaRt")
mapping <- c(intact$V1, intact$V2)
mapping <- unique(mapping)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
uniprot <- getBM(attributes= c("hgnc_symbol", "uniprotswissprot"),
                 values = mapping,
                 mart = ensembl)
names(uniprot) <- c("Prot1", "V1" )
intact <- merge(intact, uniprot, by="V1")
names(uniprot) <- c("Prot2", "V2" )
intact <- merge(intact, uniprot, by="V2")
intact <- dplyr::select(intact, Prot1, Prot2)
intact <- na.omit(intact)
intact$Source <- "intact"
write.csv(intact, "../5_PPI_API/cleaned/intact.csv")

bio1 <- read.delim("../5_PPI_API/DN/BioPlex/BioPlex_293T_Network_10K_Dec_2019.tsv")
bio2 <- read.delim("../5_PPI_API/DN/BioPlex/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv")
bio <- rbind(bio1, bio2)
rm(bio1, bio2)
# Retrieve a list of seeds 
peptides <- read.table("../5_PPI_API/DN/StringAPI/peptides.txt", quote="\"", comment.char="")
bio <- bio[, names(bio) %in% c("SymbolA", "SymbolB")]
bio2 <- bio[(bio$SymbolA %in% peptides$V1) | (bio$SymbolB %in% peptides$V1), ]
list <- unique(c(as.character(bio2$SymbolB), as.character(bio2$SymbolA)))
bio3 <- bio[(bio$SymbolA %in% list) & (bio$SymbolB %in% list), ]
list2 <- unique(c(as.character(bio3$SymbolB), as.character(bio3$SymbolA)))
bio <- bio3

bio$Source <- "bioplex"
write.csv(bio, "../5_PPI_API/cleaned/bioplex.csv")
rm(list = ls())
#################################################
string <- read.csv("../5_PPI_API/cleaned/string.csv", row.names=1)
intact <- read.csv("../5_PPI_API/cleaned/intact.csv", row.names=1)
IID_PPIs <- read.csv("../5_PPI_API/cleaned/IID_PPIs.csv", row.names=1)
biogrid <- read.csv("../5_PPI_API/cleaned/biogrid.csv", row.names=1)
bioplex <- read.csv("../5_PPI_API/cleaned/bioplex.csv", row.names=1)
huri <- read.csv("../5_PPI_API/cleaned/huri.csv", row.names=1)

ppi <- rbind(string, biogrid)
ppi <- rbind(ppi, IID_PPIs)
ppi <- rbind(ppi, intact)
ppi <- rbind(ppi, huri)
names(bioplex) <- names(ppi)
ppi <- rbind(ppi, bioplex)

rm(list=setdiff(ls(), "ppi"))

# Remove mixed labelling
ppi <- as.data.frame(ppi)
ppi$Prot1 <- as.character(ppi$Prot1)
ppi$Prot2 <- as.character(ppi$Prot2)
ppi$Source <- as.character(ppi$Source)
rownames(ppi) <- NULL



# All in caps
ppi$Prot1 <- toupper(ppi$Prot1)
ppi$Prot2 <- toupper(ppi$Prot2)
ppi$Source <- toupper(ppi$Source)

write.csv(ppi, "../5_PPI_API/cleaned/ppi.csv")
rm(list = ls())

ppi <- read.csv("../5_PPI_API/cleaned/ppi.csv", row.names = 1)

a <- as.character(ppi$Prot1)
a[nchar(a)==max(nchar(a))]

a <- unique(as.character(ppi$Prot2))
a[nchar(a)==max(nchar(a))]

ppi$Source <- NULL
ppi <- ppi[!duplicated(ppi),]
######################################################################################
# Convert aliases
library("AnnotationDbi")
library("DBI")
library("org.Hs.eg.db")
dbCon <- org.Hs.eg_dbconn()
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
aliasSymbol <- aliasSymbol[, c("alias_symbol", "symbol")]

# To consolidate hits with different name, map everything possible to entrez 
names(aliasSymbol) <- c("Prot1", "con1")
ppi = merge(ppi, aliasSymbol[!duplicated(aliasSymbol$Prot1),], by="Prot1" )
names(aliasSymbol) <- c("Prot2", "con2")
ppi = merge(ppi, aliasSymbol[!duplicated(aliasSymbol$Prot2),], by="Prot2" )

write.csv(ppi, "../5_PPI_API/cleaned/remappedppi.csv", quote = F, row.names = F)

#############################################################
################## Topological analyses #####################
#############################################################
ppi <- read.csv("../5_PPI_API/cleaned/remappedppi.csv")
ppi <- ppi[, 3:4]
ppi <- ppi[!duplicated(ppi),]
ppi <- na.omit(ppi)
ppi$con2 <- as.character(ppi$con2)
ppi$con1 <- as.character(ppi$con1)
length(unique(c(ppi$con1, ppi$con2)))
library("igraph")
g <- graph_from_data_frame(ppi, directed = F)
# Remove again multiple edge between the same nodes and self-loops
g <- simplify(g)
# Remove terminal
Isolated = which(degree(g)<2)
g = delete.vertices(g, Isolated)

dg <- decompose.graph(g) 
g<- dg[[1]]
is.directed(g) # Check that the network is not directed 
ecount(g)# Number of edges   155937
vcount(g) # Number of nodes   5705

adjacency_matrix <- igraph::as_adjacency_matrix(g, type="both")

cyto <- as.data.frame(as_edgelist(g, names = TRUE))  
write.csv(cyto, "../5_PPI_API/cleaned/cytoscapeedges.csv", quote = F, row.names = F)
 

#########################################################################################
# Select only the big continent. 
deg <-degree(g, mode="all", normalized=F)
betw <- betweenness(g, directed = FALSE, weights=NA, normalize=F)

deg <- as.data.frame(deg)
betw <- as.data.frame(betw)

merge_degree <- merge(deg, betw, by=0)

merge_degree <- merge_degree[order(-merge_degree$deg), ]
merge_degree$degrank <- 1:nrow(merge_degree)

merge_degree <- merge_degree[order(-merge_degree$betw), ]
merge_degree$betrank <- 1:nrow(merge_degree)



hubs <- merge_degree[merge_degree$deg> quantile(merge_degree$deg, probs = 0.95), ]
bottle <- merge_degree[merge_degree$betw> quantile(merge_degree$betw, probs = 0.95), ]
degs <- read.csv("../3results/corrected.csv", row.names=1)

degs <- as.character(degs[degs$pBH<0.05 & abs(degs$estimate)>log(1.2) ,]$Gene)
length(intersect(hubs$Row.names, bottle$Row.names))
central <- setdiff(intersect(hubs$Row.names, bottle$Row.names), degs)

central <- list( 
  centrality_measures = merge_degree, 
  hubs = hubs,
  bottlenecks = bottle,
  degs = degs,
  central = setdiff(intersect(hubs$Row.names, bottle$Row.names), degs)
  )

length(setdiff(intersect(hubs$Row.names, bottle$Row.names), degs)) # 202
length(intersect(hubs$Row.names, bottle$Row.names)) # 207
length(intersect(intersect(hubs$Row.names, bottle$Row.names), degs)) # 5

write.csv(merge_degree, "../5_PPI_API/cleaned/nodes_centrality.csv")
write.csv(merge_degree, "../5_PPI_API/cleaned/nodes_centrality.csv")
library("openxlsx")
write.xlsx(central, file = "../5_PPI_API/cleaned/ST4_nodes_centrality.xlsx", overwrite = T)
