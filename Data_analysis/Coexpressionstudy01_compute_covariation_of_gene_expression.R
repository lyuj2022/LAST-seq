library(tidyverse)

# load counts matrix (batch effect corrected)---------------------------------------
setwd("./input/")
wt <- read.csv("148adjusted.csv",row.names = 1)
wt <- wt[!grepl("^ERCC-",rownames(wt)), ]

# normalize total UMI
# make a function for normalization, set 100000 as scaling up factors.
CPM <- function (expr_mat) {
  norm_factor <- colSums(expr_mat,na.rm = T)
  return(100000*t(t(expr_mat)/norm_factor))
}

# normalize count by libray size
nwt <- CPM(wt)
# keep genes with average expression level >= 2. Here, 2 is an arbitrary number. 
# Genes with high expression level are usually less affected by technical noise.
nwt  <- nwt[rowSums(nwt)/148 >=2, ]

# make a function to compute cv
cv <-  function(x) {
  Cv = sd(x,na.rm = T) / mean(x,na.rm =T)
  return(Cv)
}

# compute mean and cv in case they are needed somewhere. 
#In fact, We didn't need mean and cv in the downstream analysis.
countCV <- apply(nwt, MARGIN = 1, cv)
countmean <- apply(nwt, MARGIN = 1, mean)
countCV <- unlist(countCV)
countmean <- unlist(countmean)
mean_cv_wt <-as.data.frame(cbind(countmean,countCV)) 

# load annotation including genename 
genename <- read.csv("lastsp.gene_names.txt",sep = "\t")
mean_cv_wt$gene_id <- rownames(mean_cv_wt)

# annotate genes with their names
mean_cv_wt_f <- inner_join(mean_cv_wt,genename,by="gene_id")

# load genes with TAD coordinates
TADs <- read.delim("wt_gene_TADs.bed",header = F)
colnames(TADs) <- c("Chr","st","ed","gene_name","length","tads")

# select genes located in TADs
mean_cv_TDAs <- inner_join(mean_cv_wt_f,TADs, by="gene_name")
mean_cv_TDAs <- mean_cv_TDAs %>%
  filter(tads != "")
tads_stat <- as.data.frame(table(mean_cv_TDAs$tads))

# select TADs containing >= 2 genes
tadselec <- tads_stat %>%
  filter(Freq >=2)

# retrieve TADs ID
selec_TADs <- tadselec$Var1
#retrieve gene ID
mean_cv_tads_selec <- mean_cv_TDAs[mean_cv_TDAs$tads %in% selec_TADs, ]

# retrieve counts for selected genes
nwt_ <- as.data.frame(nwt)
nwt_$gene_id <- rownames(nwt_)
clean_TAD_counts <- inner_join(mean_cv_tads_selec,nwt_,by="gene_id")
clean_TAD_counts <- clean_TAD_counts[, c(4,9:157)]

# make a function to combine SCC and p-value

flat_cor_mat <- function(cor_r, cor_p) {
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_r <- cor_r %>% filter(cor != "NA")
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# compute SCCs of gene pairs for each TADs
inde <- as.character(tadselec$Var1)
SCClist <- vector("list",length = 551)


for (i in seq_along(inde)) {
  tepID <- filter(clean_TAD_counts, tads ==inde[i])
  rownames(tepID) <- tepID$gene_name
  tepID <-tepID[,-c(1,2)]
  tepIDcor <- cor(t(tepID),method = "spearman")
  tepIDcorp <-cor_pmat(t(tepID),method = "spearman")
  #move duplicates and 1
  tepIDcor[lower.tri(tepIDcor, diag = TRUE)] <- NA
  tepIDcorp[lower.tri(tepIDcorp, diag = TRUE)] <- NA
  
  SCClist[[i]] <- flat_cor_mat(tepIDcor,tepIDcorp)
}

# merge all elements in list by row
SCC_p <- do.call(rbind, SCClist)

# the same strategy applies to TAD boundaries--------------------------------------------
TADbound <- read.delim("wt_gene_TADbound_0.05.bed",header = F)
colnames(TADbound) <- c("Chr","st","ed","gene_name","length","tads")
mean_cv_TDAbound <- inner_join(mean_cv_wt_f,TADbound, by="gene_name")

mean_cv_TDAbound <- mean_cv_TDAbound %>%
  filter(tads != "")
tads_statbound<- as.data.frame(table(mean_cv_TDAbound$tads))
tadboundselec <- tads_statbound %>%
  filter(Freq >=2)
selec_TADbound <- tadboundselec$Var1

mean_cv_tadbound_selec <- mean_cv_TDAbound[mean_cv_TDAbound$tads %in% selec_TADbound, ]

clean_TADbound_counts <- inner_join(mean_cv_tadbound_selec,nwt_,by="gene_id")
clean_TADbound_counts <- clean_TADbound_counts[, c(4,9:157)]

inde <- as.character(tadboundselec$Var1)
SCClistbound <- vector("list",length = 249)

for (i in seq_along(inde)) {
  tepID <- filter(clean_TADbound_counts, tads ==inde[i])
  rownames(tepID) <- tepID$gene_name
  tepID <-tepID[,-c(1,2)]
  tepIDcor <- cor(t(tepID),method = "spearman")
  tepIDcorp <-cor_pmat(t(tepID),method = "spearman")
  #move duplicates and 1
  tepIDcor[lower.tri(tepIDcor, diag = TRUE)] <- NA
  tepIDcorp[lower.tri(tepIDcorp, diag = TRUE)] <- NA
  SCClistbound[[i]] <- flat_cor_mat(tepIDcor,tepIDcorp)
}

SCC_pbound <- do.call(rbind, SCClistbound)

## compute SCC for all genes------------------------------------------------
allgenecounts <- inner_join(mean_cv_wt_f,nwt_,by="gene_id")
tepID <- allgenecounts[, 4:152]

rownames(tepID) <- tepID$gene_name
tepID <-tepID[, -1]

tepIDcor <- cor(t(tepID),method = "spearman")
tepIDcorp <-cor_pmat(t(tepID),method = "spearman")
# move duplicates and 1
tepIDcor[lower.tri(tepIDcor, diag = TRUE)] <- NA
tepIDcorp[lower.tri(tepIDcorp, diag = TRUE)] <- NA


cor_r <- rownames_to_column(as.data.frame(tepIDcor), var = "row")
cor_r <- gather(cor_r, column, cor, -1)
cor_r <- cor_r %>%
  filter(cor != "NA")

cor_p <- rownames_to_column(as.data.frame(tepIDcorp), var = "row")
cor_p <- gather(cor_p, column, p, -1)
cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
