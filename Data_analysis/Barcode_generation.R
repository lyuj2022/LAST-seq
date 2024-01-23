library("DNABarcodes")

# generate 6nt barcodes with 3 hamming distances ------------------------------------------------------
mySetDist6 <- create.dnabarcodes(6, dist=3, heuristic="ashlock", filter.triplets=FALSE)
mySetDist6 <- as.data.frame(mySetDist6)
write.csv(mySetDist6, file = "barcode_1c.csv", row.names = F)

# barcodes containing GAG, AAG, GGG, GCG at the 3’end are removed.

# Analysing a Set of DNA Barcodes to calculate hamming distance 
# barcodesetV <- c("AGCTAG","AGCTTC","CATGCA","CAGATC","TCACAG","GTCTAG","GTTGCA","GTGACT","ACTCGT","TGCAGA")
# analyse.barcodes(barcodesetV)

# join barcodes to LAST-primers -------------------------------------------------------------------
library("tidyverse")

barc <- read_csv("barcode_1c.csv", col_names = F )
barc_v <- as.character(barc$X1)
cobarc <- character(length(barc_v))
indices <- 1:length(barc_v)

for (i in indices) {
  conb <- paste("GGGTGAATGAAAGGTGGCTCTAATACGACTCACTATAGGGAGACGTGTGCTCTTCCGATCTNNNNNNNN", barc_v[i], "TCGTATTAGAG")
  cobarc[i] <- conb
}

cobarcd <- as.data.frame(cobarc)
write_delim(cobarcd,"barcodingprimer80.cvs")







