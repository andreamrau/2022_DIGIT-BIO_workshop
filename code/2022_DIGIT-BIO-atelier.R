## ----setup, include=FALSE-------------------------------------
knitr::opts_chunk$set(echo = TRUE)
# knitr::purl("2022_DIGIT-BIO-atelier.Rmd")


## ---- eval = FALSE--------------------------------------------
## install.packages("BiocManager")
## BiocManager::install(c("mclust", "clValid", "fpc", "ggplot2",
##                        "coseq", "clusterProfiler", "org.Mm.eg.db"))


## ---- eval=FALSE----------------------------------------------
## install.packages("remotes")
## remotes::install_github('YuLab-SMU/ggtree')


## ---- warning=FALSE, message=FALSE----------------------------
library(clValid)
library(fpc)
library(mclust)
library(coseq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)


## ---- eval=FALSE----------------------------------------------
## for(i in c("Fu-GSE60450_dds.rds", "coexp_gmm.rds", "coexp_km.rds")) {
##   FileURL <- paste0("https://github.com/andreamrau/",
##                     "2022_DIGIT-BIO_workshop/tree/main/data/", i)
##   download.file(FileURL, i)
## }


## -------------------------------------------------------------
dds <- readRDS("Fu-GSE60450_dds.rds")
coexp_gmm <- readRDS("coexp_gmm.rds")
coexp_km <- readRDS("coexp_km.rds")
set.seed(12345)


## -------------------------------------------------------------
coexp_gmm
plot(coexp_gmm, conds = dds$Status,
     collapse_reps = "average", graph = "boxplots")
plot(coexp_gmm, graphs = "probapost_boxplots")


## -------------------------------------------------------------
coexp_km
plot(coexp_km, conds = dds$Status,
     collapse_reps = "average", graph = "boxplots")
coseq::plot(coexp_km, graphs = "probapost_boxplots")


## ---- cache=TRUE----------------------------------------------
sil_gmm <- cluster::silhouette(clusters(coexp_gmm), dist(tcounts(coexp_gmm)))
summary(sil_gmm)
sil_df_gmm <- data.frame(as(sil_gmm, "matrix"))
sil_df_gmm$Cluster <- factor(sil_df_gmm$cluster)
sil_df_gmm <- sil_df_gmm[order(sil_df_gmm$cluster, sil_df_gmm$sil_width),]
ggplot(sil_df_gmm) +
  geom_col(aes(x=1:nrow(sil_df_gmm), y=sil_width, fill = Cluster)) + 
  ylab("Silhouette width") + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## ---- cache=TRUE----------------------------------------------
sil_km <- cluster::silhouette(clusters(coexp_km), dist(tcounts(coexp_km)))
summary(sil_km)
sil_df_km <- data.frame(as(sil_km, "matrix"))
sil_df_km$Cluster <- factor(sil_df_km$cluster)
sil_df_km <- sil_df_km[order(sil_df_km$cluster, sil_df_km$sil_width),]
ggplot(sil_df_km) +
  geom_col(aes(x=1:nrow(sil_df_km), y=sil_width, fill = Cluster)) + 
  ylab("Silhouette width") + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## ---- cache=TRUE----------------------------------------------
coexp_gmm_dunn <- cluster.stats(d = dist(tcounts(coexp_gmm)), 
                             clusters(coexp_gmm))$dunn
coexp_km_dunn <- cluster.stats(d = dist(tcounts(coexp_km)), 
                                clusters(coexp_km))$dunn
cat("DI (GMM) =", coexp_gmm_dunn, "\nDI (K-means) =", coexp_km_dunn, "\n")


## ---- cache=TRUE----------------------------------------------
coexp_gmm_stability <- 
  clValid(as.matrix(tcounts(coexp_gmm)), 
          nClust=max(clusters(coexp_gmm)), maxitems = nrow(coexp_gmm),
          clMethods = "model",
          validation = "stability")
summary(coexp_gmm_stability)


## ---- cache=TRUE----------------------------------------------
coexp_km_stability <- 
  clValid(as.matrix(tcounts(coexp_km)), 
          nClust=max(clusters(coexp_km)), maxitems = nrow(coexp_km),
          clMethods = "kmeans",
          validation = "stability")
summary(coexp_km_stability)


## ---- cache=TRUE----------------------------------------------
coexp_gmm_boot <- 
  clusterboot(tcounts(coexp_gmm), B=20, k=max(clusters(coexp_gmm)), 
              seed=20,bootmethod="boot",
              clustermethod=noisemclustCBI)
print(coexp_gmm_boot)


## ---- cache=TRUE----------------------------------------------
coexp_km_boot <- 
  clusterboot(tcounts(coexp_km), B=20, k=max(clusters(coexp_km)), 
              seed=20,bootmethod="boot",
              clustermethod=kmeansCBI)
print(coexp_km_boot)


## ---- cache=TRUE----------------------------------------------
adjustedRandIndex(clusters(coexp_gmm), clusters(coexp_km))

chr_info <- rowData(dds)[match(rownames(coexp_gmm), 
                               rownames(rowData(dds))),"TXCHROM"]
table(chr_info)
adjustedRandIndex(clusters(coexp_gmm), chr_info)
adjustedRandIndex(clusters(coexp_km), chr_info)


## ---- cache=TRUE----------------------------------------------
coexp_gmm_bio <- clValid(as.matrix(tcounts(coexp_gmm)), 
                        nClust = max(clusters(coexp_gmm)),
                        maxitems = nrow(coexp_gmm),
                        clMethods = c("model"),
                        validation = c("biological"),
                        annotation = "org.Mm.eg.db", 
                        GOcategory = "BP")


## ---- cache=TRUE----------------------------------------------
coexp_km_bio <- clValid(as.matrix(tcounts(coexp_km)), 
                       nClust = max(clusters(coexp_km)),
                       maxitems = nrow(coexp_km),
                       clMethods = c("kmeans"),
                       validation = c("biological"),
                       annotation = "org.Mm.eg.db", 
                       GOcategory = "BP")


## ---- cache=TRUE----------------------------------------------
optimalScores(coexp_gmm_bio)
optimalScores(coexp_km_bio)



## ---- cache=TRUE----------------------------------------------
cl <- coexp_gmm
kchoose <- 10
ego <- enrichGO(gene = rownames(cl)[which(clusters(cl) == kchoose)],
                universe = rownames(cl),
                OrgDb = org.Mm.eg.db,
                ont = "BP")


## ---- fig.height = 10-----------------------------------------
dotplot(ego, showCategory=15)


## ---- warning=FALSE-------------------------------------------
head(summary(ego))


## ---- eval=FALSE----------------------------------------------
## library(RnaSeqGeneEdgeRQL)
## library(DESeq2)
## library(Mus.musculus)
## 
## 
## ## Read in targets file --------------------------------------------------------
## targetsFile <-
##   system.file("extdata", "targets.txt", package="RnaSeqGeneEdgeRQL")
## targets <- read.delim(targetsFile, stringsAsFactors=FALSE)
## targets$ID <- rownames(targets)
## targets$Status <- factor(targets$Status,
##                          levels = c("virgin", "pregnant", "lactating"))
## targets
## 
## 
## ## Download gene-wise read counts from NCBI ------------------------------------
## FileURL <- paste(
##      "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
##      "format=file",
##      "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
##      sep="&")
## download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")
## 
## GenewiseCounts <-
##   read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",
##              row.names="EntrezGeneID")
## colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)
## counts <- GenewiseCounts[,-1] ## First column is gene lengths
## counts <- counts[-which(rowSums(counts) == 0),]
## 
## dim(counts)
## head(counts)
## 
## 
## ## Gene names ------------------------------------------------------------------
## geneid <- rownames(counts)
## genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
##                 keytype="ENTREZID")
## ## Remove multi-mapped genes
## genes <- genes[!duplicated(genes$ENTREZID),]
## all.equal(unlist(genes$ENTREZID), rownames(counts))
## dim(genes)
## head(genes)
## 
## 
## ## Differential analysis to screen DE genes ------------------------------------
## dds <- DESeqDataSetFromMatrix(counts, targets, design = ~ Status)
## rowData(dds) <- genes
## ## Remove weakly expressed genes
## keep <- rowSums(counts(dds) >= 20) >= 2
## dds <- dds[keep,]
## ## Remove genes missing a symbol or chromosome
## keep2 <- which(rowSums(is.na(rowData(dds))) == 0)
## dds <- dds[keep2,]
## dds <- DESeq(dds, test="LRT", reduced = ~1)
## res <- results(dds, alpha = 0.05)
## hist(res$pvalue, breaks = 50, main = "Raw p-values")
## 
## 
## ## Run coseq for co-expression analysis ----------------------------------------
## coexp_gmm_all <- vector("list", length = 20)
## for(i in 1:20) {
##     coexp_gmm_all[[i]] <- coseq(dds, alpha = 0.05, K = 2:20,
##                                 model = "Normal", transformation = "arcsin")
## }
## ICL_all <- unlist(lapply(coexp_gmm_all, function(x) min(ICL(x))))
## best_gmm_result <- which.min(ICL_all)
## coexp_gmm <- coexp_gmm_all[[best_gmm_result]]
## 
## coexp_km <- coseq(dds, alpha = 0.05, K = 2:20, seed=12345,
##                   model = "kmeans", transformation = "arcsin")
## 
## ## Save DE and coseq results ---------------------------------------------------
## saveRDS(dds, file = "Fu-GSE60450_dds.rds")
## saveRDS(coexp_gmm, file = "coexp_gmm.rds")
## saveRDS(coexp_km, file = "coexp_km.rds")


## ----sessionInfo, echo=FALSE----------------------------------
sessionInfo()

