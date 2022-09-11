
library(tidyverse)

######################extract bursting frequency##################
# load point estimated value of bursting frequency
bf <- read.csv("bf.csv",header = T,row.names = 1) 

# load the lower value of CI
bflow <- read.csv("bfl.csv",header = T,row.names = 1)

#discard point estimated bf without corresponding CI
bf <- cbind(bf,bflow) 
colnames(bf) <- c("X1","X2")
bf <- bf %>%
  filter(X2 != "NA")
bf$X2 <- rownames(bf)

######################extract bursting size as above #############
bz <- read.csv("bs.csv",header = T,row.names = 1)
bslow <- read.csv("bs_low.csv",header = T,row.names = 1)

bs <- cbind(bz,bslow)
colnames(bs) <- c("X1","X2")

bs <- bs %>%
  filter(X2 != "NA")
bs$X2 <- rownames(bs)

##################################################################
#combine frequency and size
bf_bs <- inner_join(bf,bs, by="X2")
rownames(bf_bs) <- bf_bs$X2
bf_bs$X2 <- NULL
colnames(bf_bs) <- c("B_freqency","B_size")

#calculate CV of  count matrix
expression <- read.csv("148adjusted.csv",header = T, row.names = 1)
transcript <- unlist(apply(expression,MARGIN = 1,mean))
transcript <- as.data.frame(transcript)
transcript$ID <- rownames(transcript)

#combine bursting kinetics and expression CV
bf_bs$ID <- rownames(bf_bs)
bf_bs_e <- inner_join(bf_bs,transcript, by= "ID")

#save R objcts
saveRDS(bf_bs_e,file = 'bf_bs_e.rds')


