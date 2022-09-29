
#load libraries -----------------------------
library(tidyverse)
library(boot) 

#make CV function-------------------------------------------------

cv <-  function(x) {
  Cv = sd(x) / mean(x)
  return(Cv)
}

#compute CV for burst frequency and size -------------------------------------------------
#load gene with TADs information
TADbound <- read.delim("wt_gene_TADbound_0.05.bed",header = F)
colnames(TADbound) <- c("Chr","st","ed","name","length","tads")

#load frequency and size
bf_bs_e <- readRDS('bf_bs_e.rds')

#merge gene location (within TADs) and burst frequency/size
fre_size_TDAbound <- inner_join(bf_bs_e,TADbound, by="name")

#discard genes without location information
fre_size_TDAbound <- fre_size_TDAbound %>%
  filter(tads != "")
tadbound_stat <- as.data.frame(table(fre_size_TDAbound$tads))

#select TADs containing >= 4 genes. This number is chosen to get as many boundaries as possible
tadboundselec <- tadbound_stat %>%
  filter(Freq >=4)

#get the median of gene number from the selected boundaries, and use it as sampling size for the control.
median(tadboundselec$Freq) # 5 genes

#subset frequency/size matrix based on selected boundaries
selec_TADbound <- tadboundselec$Var1
length(selec_TADbound)# 37 boundaries
fre_size_tadbound_selec <- fre_size_TDAbound[fre_size_TDAbound$tads %in% selec_TADbound,]

#calculate frequency and size CV of genes within the selected boundaries
fre_size_tadbound_selec_by <- fre_size_tadbound_selec %>% group_by(tads)
SD_tadbound_fs <- fre_size_tadbound_selec_by %>% summarise(CVfreq=cv(B_freqency),CVsize=cv(B_size))
SD_tadbound_fs$tads <-NULL

#bootstrapping median_frequencyCV
boot_TADbound_FreqCV <- boot(data = SD_tadbound_fs$CVfreq,statistic = function(x,i) median(x[i]),R = 10000)
boot_TADbound_FreqCV <- boot_TADbound_FreqCV$t
type <- c(rep("TADboundary",10000))
boot_TADbound_FreqCV <-data.frame(FreqCV=boot_TADbound_FreqCV,type=type)
median(boot_TADbound_FreqCV$FreqCV) #0.4907204

# generate control--------------------------------------------------------------------------------------------
# the same logic as generating control for TADs

#subset frequency and size------------------------------------------------------
fre_size <- bf_bs_e[,1:2]

#generate empty vectors
randoms <- vector("list",length = 10000)
inde <- 1:10000
freqCVbound <- vector("numeric",length = 10000)
sizeCVbound <- vector("numeric",length = 10000)
le <- 1:37
medfreqCVbound <- vector("numeric",length = 37)
medsizeCVbound <- vector("numeric",length = 37)

#generate 37 CV numbers as those generated from 37 boundaries.
for (n in le) {
  for (i in inde){
    randoms[[i]] <- sample(1:5156,4)
    subgene <-  fre_size[randoms[[i]],]
    freqCVbound[i] <- cv(subgene$B_freqency)
    sizeCVbound[i] <- cv(subgene$B_size)
  }
  medfreqCVbound[n] <- median(freqCV)
  medsizeCVbound[n] <- median(sizeCV)
}

###freqCV_boot
boot_contrl_FreqCVbound <- boot(data = medfreqCVbound,statistic = function(x,i) median(x[i]),R = 10000)
boot_contrl_FreqCVbound <- boot_contrl_FreqCVbound$t
type <- c(rep("RANDOM",10000))
boot_contrl_FreqCVbound <-data.frame(FreqCV=boot_contrl_FreqCVbound,type=type)
median(boot_contrl_FreqCVbound$FreqCV)# 0.5035513
