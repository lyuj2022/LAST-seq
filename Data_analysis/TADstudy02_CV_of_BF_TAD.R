# load libraries -----------------------------
library(tidyverse)
library(boot)

# make CV function-------------------------------------------------

cv <-  function(x) {
  Cv = sd(x) / mean(x)
  return(Cv)
}

# compute CV for burst frequency and size -------------------------------------------------
## load genes with TADs information
TADs <- read.delim("wt_gene_TADs.bed",header = F)
colnames(TADs) <- c("Chr","st","ed","name","length","tads")

## load frequency and size
bf_bs_e <- readRDS('./input/bf_bs_e.rds')

## merge gene location (in TADs) and burst frequency/size
fre_size_TDAs <- inner_join(bf_bs_e,TADs, by="name")

## discard genes without location information
fre_size_TDAs <- fre_size_TDAs %>%
  filter(tads != "")
tads_stat <- as.data.frame(table(fre_size_TDAs$tads))

## select TADs containing >= 8 genes. This number is chosen to get as many TADs as possible
tadselec <- tads_stat %>%
  filter(Freq >=8)

## get the median of gene number from the selected TADs, and use it as sampling size for the control.
median(tadselec$Freq) # 10 genes

## subset frequency and size matrix based on selected TADs
selec_TADs <- tadselec$Var1
length(selec_TADs)# 59 TADs
fre_size_tads_selec <- fre_size_TDAs[fre_size_TDAs$tads %in% selec_TADs, ]

## calculate frequency and size CV of genes within the selected TADs
fre_size_tads_selec_by <- fre_size_tads_selec %>% group_by(tads)
SD_tads_fs <- fre_size_tads_selec_by %>% summarise(CVfreq=cv(B_freqency),CVsize=cv(B_size))
SD_tads_fs$tads <- NULL

## bootstrapping median_frequencyCV
boot_TADs_FreqCV <- boot(data = SD_tads_fs$CVfreq,statistic = function(x,i) median(x[i]),R = 10000)
boot_TADs_FreqCV <- boot_TADs_FreqCV$t
type <- c(rep("TAD",10000))
boot_TADs_FreqCV <-data.frame(FreqCV=boot_TADs_FreqCV,type=type)
median(boot_TADs_FreqCV$FreqCV) #0.4697341

# generate control--------------------------------------------------------------------------------------------
                         
## how to generate the control
### 1, randomly select 10 genes from the frequency/size matrix, then compute the CV of frequency and size.
### The gene number is determined by the median of gene number of selected TADs.

### 2, repeat 1 for 10000 times and get the median of 10000 frequency/size CVs. As a result, we get one frequency CV and one size CV from 10000 times iteration. 
### since we have more than 5000 genes in the frequency/size matrix, a random selection of 10 genes is not representative. 
### To overcome this issue, we do the sampling 10000 times. For every time, we sample 10 genes and calculate the frequency CV and size CV.
### Then, we take the median of 10000 frequency/size CVs. The median represents the frequency/size CV of a virtual TAD containing 10 genes.

### 3, we repeat 1 and 2 for 59 times so that we have 59 frequency/size CVs from 59 virtual TADs.                        
                      
# subset frequency and size------------------------------------------------------
fre_size <- bf_bs_e[, 1:2]
                         
## generate empty vectors
randoms <- vector("list",length = 10000)
inde <- 1:10000
freqCV <- vector("numeric",length = 10000)
sizeCV <- vector("numeric",length = 10000)
le <- 1:59
medfreqCV <- vector("numeric",length = 59)
medsizeCV <- vector("numeric",length = 59)

## generate 59 CV numbers as those generated from 59 TADs.
for (n in le) {
  #get representative CV 
  for (i in inde) {
    #randomly select 10 genes, assuming they're in the same TAD. then get CVs
    randoms[[i]] <- sample(1:5156,10)
    subgene <-  fre_size[randoms[[i]], ]
    freqCV[i] <- cv(subgene$B_freqency)
    sizeCV[i] <- cv(subgene$B_size)
  }
  medfreqCV[n] <- median(freqCV)
  medsizeCV[n] <- median(sizeCV)
}


## freqCV_boot
boot_contrl_FreqCV <- boot(data = medfreqCV,statistic = function(x,i) median(x[i]),R = 10000)
boot_contrl_FreqCV <- boot_contrl_FreqCV$t
type <- c(rep("RANDOM",10000))
boot_contrl_FreqCV <-data.frame(FreqCV=boot_contrl_FreqCV,type=type)
median(boot_contrl_FreqCV$FreqCV)#0.5042821
