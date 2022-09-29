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
##############################################
# the principle to extend the TADs  boundary by 5% of TADs
# 1, focus on the right coordinates of TADs since they are boundaries, except the last TAD.
# 2, for each TAD, we compute the distance of 5% of the TAD. i.e. frn
# 3, for the i th boundary with location number ED, the starting coordinate is ED - frn of the i th TAD. The ending coordinate is ED + frn of the i+1 th TAD

#retrieve chromosome-separated data
flw1 <- flw[[1]]
colnames(flw1) <- c("ID","ST","ED","TADs")

#compute TAD boundary coordinates with 5% extension.
#compute the distance of 5% of TADs
flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)

#extend the left coordinate of the boundary by subtracting the the distance of 5% of TAD
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

#discard the first elements of the frn vector and fill the last elements with 0. As a resulte the i th was replaced by the i+1 th
rear <- flw1$frn
rear <-rear[-1]
rear[176] <- 0

#extend the right  coordinate of the boundary by adding  the distance of 5% of the next TAD
flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr1 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[2]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[197] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr2 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)

################################
flw1 <- flw[[3]]
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[163] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr3 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)

################################
flw1 <- flw[[4]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr4 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)

################################
flw1 <- flw[[5]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr5 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[6]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr6 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[7]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr7 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[8]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr8 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[9]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr9 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[10]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr10 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[11]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr11 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[12]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr12 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[13]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr13 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[14]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr14 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[15]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr15 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[16]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr16 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[17]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr17 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[18]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr18 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[19]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr19 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[20]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr20 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[21]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr21 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[22]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr22 <- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)
################################
flw1 <- flw[[23]]
n <- nrow(flw1)
colnames(flw1) <- c("ID","ST","ED","TADs")

flw1 <- flw1 %>%
  mutate(frn = (ED-ST)/20)
flw1 <- flw1 %>%
  mutate(nST = ED-frn)

rear <- flw1$frn
rear <-rear[-1]
rear[n] <- 0

flw1 <- cbind(flw1,rear)
flw1 <- flw1 %>%
  mutate(nED = ED+rear)
Chr23<- data.frame(flw1$ID,flw1$nST,flw1$nED,flw1$TADs)

#######meg
TADbound <- rbind(Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23)
options(scipen=999)
write.table(TADbound,"TADbound_0.05.bed",sep = '\t',row.names= F,col.names=F,quote = FALSE)
