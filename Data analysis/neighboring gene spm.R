library(tidyverse)

##load dataset to get expression matrix, gene coordinates, and gene orientations#########################
#count matrix, nwt for wt, nscc for SCC4-/-
adjusted <- read.csv("nwt.csv")

#gene coordinates
genecor <- read.csv("genecoordinate.csv")
genecor <- genecor[,-5]
colnames(genecor) <- c("CHR","ST","ED","name")

#gene orientations
gene_dir <- read.csv("gene_dirextion.csv")

#gene names, this information is not required for downsteam analysis
gene_name <- read.csv("gene_names.txt",sep = "\t")
colnames(gene_name) <- c("ID","name")

#merge all gene information
gene_dir_name <- inner_join(gene_dir,gene_name,by="ID")
gene_dir_name_coor_switch <- inner_join(genecor,gene_dir_name,by="name")

#merge count matrix and gene_dir_name_coor_switch
counts <- adjusted
counts$ID <- counts$X
counts$X <- NULL
countscor <- inner_join(gene_dir_name_coor_switch,counts)
countscor$V1 <- NULL

#split chromosomes
table(countscor$CHR)
Chr1 <- countscor %>%
  filter(CHR=="chr1")
Chr2 <- countscor %>%
  filter(CHR=="chr2")
Chr3 <- countscor %>%
  filter(CHR=="chr3")
Chr4 <- countscor %>%
  filter(CHR=="chr4")
Chr5 <- countscor %>%
  filter(CHR=="chr5")
Chr6 <- countscor %>%
  filter(CHR=="chr6")
Chr7 <- countscor %>%
  filter(CHR=="chr7")
Chr8 <- countscor %>%
  filter(CHR=="chr8")
Chr9 <- countscor %>%
  filter(CHR=="chr9")
Chr10 <- countscor %>%
  filter(CHR=="chr10")
Chr11 <- countscor %>%
  filter(CHR=="chr11")
Chr12 <- countscor %>%
  filter(CHR=="chr12")
Chr13 <- countscor %>%
  filter(CHR=="chr13")
Chr14 <- countscor %>%
  filter(CHR=="chr14")
Chr15 <- countscor %>%
  filter(CHR=="chr15")
Chr16 <- countscor %>%
  filter(CHR=="chr16")
Chr17 <- countscor %>%
  filter(CHR=="chr17")
Chr18 <- countscor %>%
  filter(CHR=="chr18")
Chr19 <- countscor %>%
  filter(CHR=="chr19")
Chr20 <- countscor %>%
  filter(CHR=="chr20")
Chr21 <- countscor %>%
  filter(CHR=="chr21")
Chr22 <- countscor %>%
  filter(CHR=="chr22")
ChrX <- countscor %>%
  filter(CHR=="chrX")

##################################calculate distance and correlation coeffient for neighboring gene pairs.####################################################################
#label gene direction. '1' represents '+'; '0' represents '-'
Chrdir <- Chr1 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  )) 

#sort gene coordintes in ascending order
Chrdir <- Chrdir %>% arrange(ST) 

#retrieve gene orientations
direction <- as.vector(Chrdir$direct)

#retrieve gene coordinates
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)

#retrieve count matrix
chr1count <-as.matrix(Chrdir[7:154]) 

#creat empty vectors
l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

#compute orientation between neighboring genes , distance and spearman's r
for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  stat = neigbstate[i]
  if (stat == 0){
    dis3[i]=ST[i+1]-ED[i]
  } else if (stat==1){
    dis3[i]=ST[i+1]-ED[i]
  } else {
    dis3[i]=ST[i+1]-ED[i]
  }
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

#merge orientation, distance and spearman's r
chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))

#define gene orientation
chr1neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))

########## repeat the same code for chromosome2-22 and X
Chrdir <- Chr2 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chr1dir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr2neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr3 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr3neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr4 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr4neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr5 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr5neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr6 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr6neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))

############################
Chrdir <- Chr7 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr7neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr8 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr8neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr9 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr9neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr10 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr10neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr11 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 


l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr11neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr12 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr12neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr13 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr13neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr14 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr14neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr15 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr15neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))

############################
Chrdir <- Chr16 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr16neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr17 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr17neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr18 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr18neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr19 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr19neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr20 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr20neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))
############################
Chrdir <- Chr21 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr21neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))

############################
Chrdir <- Chr22 %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chr22neigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))

############################
Chrdir <- ChrX %>%
  mutate(direct = case_when(
    V7=="+"~1,
    TRUE~0
  ))

Chrdir <- Chrdir %>% arrange(ST)

direction <- as.vector(Chrdir$direct)
ST <- as.numeric(Chrdir$ST)
ED <- as.numeric(Chrdir$ED)
chr1count <-as.matrix(Chrdir[7:154]) 

l <- nrow(chr1count)-1
neigbstate <- vector("numeric",length = l)
dis3 <- vector("numeric",length = l)
neigPCC <- vector("numeric",length = l)

for (i in 1:l) {
  neigbstate[i]= direction[i+1]-direction[i]
  dis3[i]=ST[i+1]-ED[i]
  neigPCC[i]=cor(chr1count[i,],chr1count[i+1,],method = "spearman")
}

chrneigborgene <- as.data.frame(cbind(neigbstate,dis3,neigPCC))
chrXneigborgene <- chrneigborgene %>%
  mutate(type=case_when(
    neigbstate==1~"divergent",
    neigbstate==0~"same",
    TRUE~"convergent"
  ))


####merge all chromosomes
allmerg <- rbind(chr1neigborgene,chr2neigborgene,chr3neigborgene,chr4neigborgene,chr5neigborgene,chr6neigborgene,chr7neigborgene,chr8neigborgene,chr9neigborgene,chr10neigborgene,chr11neigborgene,
                 chr12neigborgene,chr13neigborgene,chr14neigborgene,chr15neigborgene,chr16neigborgene,chr17neigborgene,chr18neigborgene,chr19neigborgene,chr20neigborgene,chr21neigborgene,
                 chr22neigborgene,chrXneigborgene)

#when genes are overlapping, let the distance equal '0'
allmerg <- allmerg%>%
  mutate(distance= case_when(
    dis3 < 0 ~ 0,
    TRUE~1
  ))
allmerg$distance = allmerg$dis3*allmerg$distance

#order distance
allmerg <- allmerg %>%arrange(distance)

#classify genes by their relative orientations
allconvergent <- allmerg %>%
  filter(type=='convergent')
alldivergent <- allmerg %>%
  filter(type=='divergent')
allsame <- allmerg %>%
  filter(type=='same')

##################################################### plot #####################################################
library(ggpubr)

allmergF <- allmerg%>%
  filter(distance<=1000)

table(allmergF$type)

my_comparisons <- list(c('convergent','divergent'),c('divergent','same'))
p <- allmergF%>% ggplot(aes(x=type,y=neigPCC,color=type))+
  geom_boxplot(lwd=0.25,width=0.75, fatten=1,outlier.size = 0.1)+
  stat_compare_means(comparisons =my_comparisons,size=1.5)+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.position='none',
    legend.title=element_blank(),
    legend.key=element_blank(),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25))+
  xlab("")+
  ylab("Spearman r")
ggsave(filename = "neighboringgeneSPM_1000.pdf",plot = p,width=65,height = 45,units = "mm", dpi = 300)
################################

allmergF <- allmerg%>%
  filter(distance<=10000)

table(allmergF$type)

my_comparisons <- list(c('convergent','divergent'),c('divergent','same'))
p <- allmergF%>% ggplot(aes(x=type,y=neigPCC,color=type))+
  geom_boxplot(lwd=0.25,width=0.75, fatten=1,outlier.size = 0.1)+
  stat_compare_means(comparisons =my_comparisons,size=1.5)+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.position='none',
    legend.title=element_blank(),
    legend.key=element_blank(),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25))+
  xlab("")+
  ylab("Spearman r")
ggsave(filename = "neighboringgeneSPM_10000.pdf",plot = p,width=65,height = 45,units = "mm", dpi = 300)
########################################################################
allmergF <- allmerg%>%
  filter(distance<=50000)

table(allmergF$type)


p <- allmergF%>% ggplot(aes(x=type,y=neigPCC,color=type))+
  geom_boxplot(lwd=0.25,width=0.75, fatten=1,outlier.size = 0.1)+
  stat_compare_means(comparisons =my_comparisons,size=1.5)+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.position='none',
    legend.title=element_blank(),
    legend.key=element_blank(),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25))+
  xlab("")+
  ylab("Spearman r")
ggsave(filename = "neighboringgeneSPM_50K.pdf",plot = p,width=65,height = 45,units = "mm", dpi = 300)
