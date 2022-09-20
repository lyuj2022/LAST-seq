library(tidyverse)

#load TADs coordinates 
TAD <- read.csv("wt_min300k_max1M_step100k_thres0.05_delta0.01_fdr_domains.bed",sep = "\t",header = F)
TADcoor <- TAD[1:4]

#split data by Chromosome N and deposit them in a list
Group <- unique(TADcoor$V1)
ind <- 1:length(Group)

flw <- vector("list", length(Group))

for (i in ind ) {
  dat <- TADcoor%>%
    filter(V1==Group[i])
  flw[[i]] <- dat
}

#######################################################
# to minimize the effect of TAD boundaries, the distance of 5% of TAD  was subtracted from the TAD 
# 1, compute the distance of 5% of TAD, i.e. frn
# 2, for each TAD, the new ending coordinate is the right coordinate minus frn. The new starting coordinate is the new ending coordinate minus 18 times frn.

flw1 <- flw[[1]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr1b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[2]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr2b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[3]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr3b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[4]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr4b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[5]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr5b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[6]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr6b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[7]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr7b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[8]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr8b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[9]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr9b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[10]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr10b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[11]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr11b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[12]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr12b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[13]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr13b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[14]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr14b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[15]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr15b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[16]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr16b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[17]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr17b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[18]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr18b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[19]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr19b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[20]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr20b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[21]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr21b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[22]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr22b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)
########################
flw1 <- flw[[23]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(bED = ED-frn)
flw1 <- flw1 %>%
  mutate(bST = bED-18*frn)
Chr23b <- data.frame(flw1$ID,flw1$bST,flw1$bED,flw1$TADs)

TADmiddle <- rbind(Chr1b,Chr2b,Chr3b,Chr4b,Chr5b,Chr6b,Chr7b,Chr8b,Chr9b,Chr10b,Chr11b,Chr12b,Chr13b,Chr14b,Chr15b,Chr16b,Chr17b,Chr18b,Chr19b,Chr20b,Chr21b,Chr22b,Chr23b)
options(scipen=999)
write.table(TADmiddle,"TADmiddle_0.05.bed",sep = '\t',row.names= F,col.names=F,quote = FALSE)
