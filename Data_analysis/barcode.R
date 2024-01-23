library("DNABarcodes")

#generate barcode with 6nt,3 hamming distance, with the Illumina's Sequencing by Synthesis technology, the triplet filter is unnecessary
mySetDist6 <- create.dnabarcodes(6, dist=3, heuristic="ashlock", filter.triplets=FALSE)
mySetDist6 <- as.data.frame(mySetDist6)
write.csv(mySetDist6, file = "barcode_1c.csv", row.names = F)

#6nt, 2 hamming distance
mySetDist6_2 <- create.dnabarcodes(6, dist=2, heuristic="ashlock", filter.triplets=FALSE)
mySetDist6_2 <- as.data.frame(mySetDist6)
write.csv(mySetDist6_2, file = "barcode_2hamming.csv", row.names = F)


#Analysing a Set of DNA Barcodes
barcodesetV <- c("AGCTAG","AGCTTC","CATGCA","CAGATC","TCACAG","GTCTAG","GTTGCA","GTGACT","ACTCGT","TGCAGA")
analyse.barcodes(barcodesetV)

tmp <- as.character(mySetDist6$mySetDist6)
analyse.barcodes(tmp)

#################generate barcoding primers #####################################
library("tidyverse")

barc <- read_csv("barcode_1c.csv", col_names = F )
cobarc <- character(length(barc_v))
barc_v <- as.character(barc$X1)
indices <- 1:length(barc_v)

for (i in indices) {
  conb <- paste("GGGTGAATGAAAGGTGGCTCTAATACGACTCACTATAGGGAGACGTGTGCTCTTCCGATCTNNNNNNNN", barc_v[i], "TCGTATTAGAG")
  cobarc[i] <- conb
}

cobarcd <- as.data.frame(cobarc)
write_delim(cobarcd,"barcodingprimer80.cvs")







