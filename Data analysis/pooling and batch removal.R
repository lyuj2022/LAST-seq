#rm(list=ls()) 
library(tidyverse)

#load the matrix containing read mapping info and rename the cell ID for plate1
pool1read <- read.csv("pool1.readspercell.txt",sep = "\t",header = T)
pool1read$RG <-  paste(1,pool1read$RG,sep = "_")

pool1read <- pool1read %>% pivot_wider(names_from = RG, values_from = N)
pool1read$'1_bad' <- NULL

#load the matrix containing read mapping info and rename the cell ID for plate2
pool2read <- read.csv("pool2.readspercell.txt",sep = "\t",header = T)
pool2read$RG <-  paste(2,pool2read$RG,sep = "_")

pool2read <- pool2read %>% pivot_wider(names_from = RG, values_from = N)
pool2read$'2_bad' <- NULL

#merge pool1 and pool2
poolread_1_2 <- inner_join(pool1read,pool2read,by="type")
poolread_1_2 <- as.data.frame(poolread_1_2)
rownames(poolread_1_2) <- poolread_1_2$type
poolread_1_2$type <- NULL

poolread_1_2_T <-as.data.frame(t(poolread_1_2)) 
poolread_1_2_T$total <- rowSums(poolread_1_2_T)

#compute mapping rate of exon 
poolread_1_2_T <- poolread_1_2_T %>%
  mutate(exonrate = 100*Exon/total)

#keep cells with exon% > 50%
poolread_1_2_f <- poolread_1_2_T %>%
  filter(exonrate > 50)

#keep cells with total reads >=100000
poolread_1_2_f <- poolread_1_2_f %>%
  filter(total >= 100000) 

#get the cell index
reads_f <- rownames(poolread_1_2_f)

#######################plot exon reads distribution#########################################
library(scales)
poolexonreads <- as.data.frame(poolread_1_2_f$Exon)

#make a function to change scientific format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

colnames(poolexonreads) <- 'Exonreads'
p <- poolexonreads %>% ggplot(aes(x=" ",y=Exonreads))+
  geom_jitter(size=0.,width = 0.25)+
  geom_boxplot(lwd=0.25, outlier.shape = NA, width=0.5, fatten=1,alpha=0,color='red')+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.title=element_blank(),
    legend.key=element_blank(),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25),
  )+
  scale_y_log10(limits=c(1e5, 1e8), label=scientific_10) +
  annotation_logticks(sides = "bl", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))+
  xlab("eHAP cells")+
  ylab("reads mapped to exons")

ggsave(filename = "reads distribution.pdf",plot = p,width=45/2,height = 45,units = "mm", dpi = 300)
###########################################################################################################################

###########################################################################################################################
#load count matrix for plate1 and plate2
pool1 <- readRDS("pool1.dgecounts.rds")
pool2 <- readRDS("pool2.dgecounts.rds")

#use the UMI counts, exclude reads mapped to ERCC and renane the cell ID
pool1_UMI <- pool1$umicount$exon$all
pool1_UMI_ex <- as.matrix(pool1_UMI)
pool1_UMI_ex_extr <- pool1_UMI_ex[!grepl("^ERCC-",rownames(pool1_UMI_ex )),]
colnames(pool1_UMI_ex_extr) <- paste(1,colnames(pool1_UMI_ex_extr),sep = "_")

pool2_UMI <- pool2$umicount$exon$all
pool2_UMI_ex <- as.matrix(pool2_UMI)
pool2_UMI_ex_extr <- pool2_UMI_ex[!grepl("^ERCC-",rownames(pool2_UMI_ex )),]
colnames(pool2_UMI_ex_extr) <- paste(2,colnames(pool2_UMI_ex_extr),sep = "_")

pool1_UMI_ex_extr <- as.data.frame(pool1_UMI_ex_extr)
pool1_UMI_ex_extr$ID  <- rownames(pool1_UMI_ex_extr)

pool2_UMI_ex_extr <- as.data.frame(pool2_UMI_ex_extr)
pool2_UMI_ex_extr$ID  <- rownames(pool2_UMI_ex_extr)

#merge plate1 and plate2
mat <- inner_join(pool1_UMI_ex_extr,pool2_UMI_ex_extr)

#subset the selected cells
rownames(mat) <- mat$ID
mat$ID <- NULL
mat[is.na(mat)] <- 0
mat <- mat[,colnames(mat) %in% reads_f]
##################################################################################################

##################################################################################################
#check the batch effect by PCA
megM <- as.matrix(mat)
M <- t(megM)

#remove constant/zero column
Mc <- M[ , which(apply(M, 2, var) != 0)]

##PCA
pca <- prcomp(Mc, scale = TRUE)

#plot by PC1 and PC2
plot(pca$x[,1],pca$x[,2])

#the percentage of variance accounted by PC

pca.var <-  pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="screen plot",xlab = "PC",ylab = "percent variation")

#reorganise data and plot by ggplot2
pca.data <- data.frame(sample=rownames(pca$x),X=pca$x[,1],Y=pca$x[,2])
pca.data$label <- c(rep("pool1",73),rep("pool2",75))

p <- ggplot(data = pca.data, aes(x=X,y=Y,color= label))+
  geom_point(size=.1)+
  xlab(paste("PC1-",pca.var.per[1], "%", sep = ""))+
  ylab(paste("PC2-",pca.var.per[2], "%", sep = ""))+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.title=element_blank(),
    legend.key=element_blank(),
    legend.text = element_text(size=6),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25)
    
  )

ggsave(filename = "before batch correction.pdf",plot = p,width=45*2,height = 45*2,units = "mm", dpi = 300)
#######################################################################################################################

#######################################################################################################################
# remove batch effect
library(sva)

count_matrix <- mat
colnames(count_matrix)

batch <- c(rep(1, 73), rep(2, 75))

adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)

write.csv(adjusted,"148adjusted.csv")#pool1&2
###########################################################################

###############check batch effect by PCA#################################################

megM <- as.matrix(adjusted)
M <- t(megM)

#remove constant/zero column
Mc <- M[ , which(apply(M, 2, var) != 0)]

##PCA
pca <- prcomp(Mc, scale = TRUE)

#plot using PC1 and PC2
plot(pca$x[,1],pca$x[,2])

#the percentage of variance accounted by PC

pca.var <-  pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="screen plot",xlab = "PC",ylab = "percent variation")

#reorganise data and plot by ggplot2
pca.data <- data.frame(sample=rownames(pca$x),X=pca$x[,1],Y=pca$x[,2])
pca.data$label <- c(rep("pool1",73),rep("pool2",75))
p <- ggplot(data = pca.data, aes(x=X,y=Y,color= label))+
  geom_point(size=0.1)+
  xlab(paste("PC1-",pca.var.per[1], "%", sep = ""))+
  ylab(paste("PC2-",pca.var.per[2], "%", sep = ""))+
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.title=element_blank(),
    legend.key=element_blank(),
    legend.text = element_text(size=6),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25)
    
  )

ggsave(filename = "after batch correction.pdf",plot = p,width=45*2,height = 45*2,units = "mm", dpi = 300)



