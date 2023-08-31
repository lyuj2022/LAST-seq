
#BiocManager::install("BASiCS")
library(BASiCS)
library(SingleCellExperiment)
library(tidyverse)

####################WT

wt <- read.csv("WTadjusted.csv",row.names = 1)
wt$ID <- rownames(wt)
scc4 <- read.csv("SCCadjusted.csv",row.names = 1)
scc4$ID <- rownames(scc4)
meg <- inner_join(wt,scc4,by="ID")


#meg <- meg[rowSums(meg>0)>30,]


wt <- meg[,1:149]
rownames(wt) <- wt$ID
wt$ID <- NULL
scc4 <- meg[,149:302]
rownames(scc4) <- scc4$ID
scc4$ID <- NULL

colnames(wt)
#
wt <-as.matrix(wt)
SCEwt<- SingleCellExperiment(
  assays = list(counts = wt),
  colData = data.frame(BatchInfo =c(rep(1,73),rep(2,75))))

counts(SCEwt)

#Typically, setting N=20000, Thin=20 and Burn=10000 leads to stable results.
#Please ensure the acceptance rates displayed in the console output of BASiCS_MCMC are around 0.44.
#If they are too far from this value, you should increase N and Burn
Chainwt <- BASiCS_MCMC(
  Data = SCEwt,
  N = 20000, Thin = 20, Burn = 10000,
  PrintProgress = TRUE, Regression = TRUE,WithSpikes = FALSE
)

#################SCC4
colnames(scc4)
#scc4 <- scc4[rowSums(scc4>0)>0,]
scc4 <-as.matrix(scc4)
SCEscc4<- SingleCellExperiment(
  assays = list(counts = scc4),
  colData = data.frame(BatchInfo =c(rep(1, 77), rep(2, 76))))

#Typically, setting N=20000, Thin=20 and Burn=10000 leads to stable results.
#Please ensure the acceptance rates displayed in the console output of BASiCS_MCMC are around 0.44.
#If they are too far from this value, you should increase N and Burn
Chainscc4 <- BASiCS_MCMC(
  Data = SCEscc4,
  N = 20000, Thin = 20, Burn = 10000,
  PrintProgress = TRUE, Regression = TRUE,WithSpikes = FALSE
)

##test differnetial expressed genes
scc_wt_Test <- BASiCS_TestDE(Chain1 = Chainscc4, Chain2 = Chainwt,
                      GroupLabel1 = "scc4", GroupLabel2 = "wt",
                      EpsilonM = log2(2), EpsilonD = log2(2),
                      EFDR_M = 0.05, EFDR_D = 0.05,
                      Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)


-------------------------------------------------------------
  1675 genes with a change in mean expression:
  - Higher expression in scc4 samples: 967
- Higher expression in wt samples: 708
- Fold change tolerance = 200% 
- Probability threshold = 0.80675
- EFDR = 4.99% 
- EFNR = 12.33% 
-------------------------------------------------------------
  
  -------------------------------------------------------------
  1300 genes with a change in over dispersion:
  - Higher dispersion in scc4 samples: 639
- Higher dispersion in wt samples: 661
- Fold change tolerance = 200% 
- Probability threshold = 0.89075
- EFDR = 4.99% 
- EFNR = 47.82% 
NOTE: differential dispersion assessment only applied to the 
16716 genes for which the mean did not change. 
and that were included for testing. 
--------------------------------------------------------------
  -------------------------------------------------------------
  909 genes with a change in residual over dispersion:
  - Higher residual dispersion in scc4 samples: 540
- Higher residual dispersion in wt samples: 369
- Distance tolerance = 0.41
- Probability threshold = 0.89275
- EFDR = 5% 
- EFNR = 64.28% 
NOTE: differential residual dispersion assessment applied to 
18279 genes expressed in at least 2 cells per condition 
and that were included for testing. 
#--------------------------------------------------------------

BASiCS_PlotDE(scc_wt_Test, Plots = "MA", Parameters = "Mean")





#--------

df <-scc_wt_Test$TableMean


#write.csv(df,"scc_df.csv")


colnames(df)

dff <- df %>%
  mutate(gene_type = case_when(ResultDiffMean=="wt+" ~ "down",
                               ResultDiffMean=="scc4+" ~ "up",
                               TRUE ~ "ns"))
dff <- dff %>%
  filter(Mean1>0.1 |  Mean2>0.1)

#Odf <- df %>%
  #filter(Mean1<=1|Mean2<=1)


DEG <- dff %>% 
  filter(gene_type != "ns")

write.csv(DEG,"DEG.csv")

as.numeric(nrow(dff))
table(dff$gene_type)
down <- 100*708 /16896
up <- 100*967/16896

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

saveRDS(dff,file = "SCC4_vs_WT_DGE.rds")

ggplot(data=dff,aes(x = MeanLog2FC,
                               y = ProbDiffMean,
                               color = gene_type)) + 
  geom_point(size=0.5)+
  #geom_hline(yintercept = 0.9,
             #linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1/2), log2(2)),
             linetype = "dashed") +
  scale_colour_manual(values = cols)+ # Modify point colour
  xlim(-6, 6)+
  theme_bw()+
  #remove grid and change border thickness and change tick label size and legend front size
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size=0.5),
        axis.line.y = element_line(colour = "black", size=0.5),
        axis.text.x = element_text(face="bold", color="black", 
                                   size=10),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=10),
        legend.text = element_text(size=10))+ 
  #change x and y labels
  #xlab("log10(bursting frequency p-value)")+  
  #ylab("log2(fold_change)")+ 
  #remove lengend title
  labs(color = " ")+
  #change legend size
  guides(colour = guide_legend(override.aes = list(size=1)))+
  annotate("text", x = -4, y = 0.8,
           label = "4.2%", color = "firebrick")+
  guides(colour = guide_legend(override.aes = list(size=1)))+
  annotate("text", x = 4, y = 0.8,
           label = "5.7%", color = "firebrick")





############check dispersion ###################not necessary

disper <- scc_wt_Test$TableDisp
disper <- disper[disper$GeneName %in% dff$GeneName,]

#disper %>% ggplot(aes(x=log2(Disp1),y=log2(Disp2)))+
  #geom_point(size=0.2)+
  #geom_abline(slope = 1,intercept = 0,color="red")

disperexclude <- disper %>%
  filter(ResultDiffDisp == "ExcludedFromTesting" )

1675

#disper <- disper %>%
  #filter(ResultDiffDisp != "ExcludedFromTesting" )




disperF <- disper %>%
  mutate(gene_type = case_when(DispFC >= 2 & ProbDiffDisp > 0.89 ~ "up",
                               DispFC <= 0.5 & ProbDiffDisp > 0.89 ~ "down",
                               TRUE ~ "ns"))
saveRDS(disperF,file = 'SCC4_vs_WT_dispersion.rds')



DEGDisp <- disperF %>% 
  filter(gene_type != "ns")

write.csv(DEGDisp,"DEGDisp.csv")

as.numeric(nrow(disper))
table(disperF$gene_type)
down <- 100* 686/16896
up <- 100*717/16896

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")


ggplot(data=disperF,aes(x = DispLog2FC,
                    y = ProbDiffDisp,
                    color = gene_type)) + 
  geom_point(size=0.5)+
  #@geom_hline(yintercept = 0.85,
             #linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1/2), log2(2)),
             linetype = "dashed") +
  scale_colour_manual(values = cols)+ # Modify point colour
  xlim(-10, 10)+
  theme_bw()+
  #remove grid and change border thickness and change tick label size and legend front size
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size=0.5),
        axis.line.y = element_line(colour = "black", size=0.5),
        axis.text.x = element_text(face="bold", color="black", 
                                   size=10),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=10),
        legend.text = element_text(size=10))+ 
  #change x and y labels
  #xlab("log10(bursting frequency p-value)")+  
  #ylab("log2(fold_change)")+ 
  #remove lengend title
  labs(color = " ")+
  #change legend size
  guides(colour = guide_legend(override.aes = list(size=1)))+
  annotate("text", x = -8, y = 0.8,
           label = "4.1%", color = "firebrick")+
  guides(colour = guide_legend(override.aes = list(size=1)))+
  annotate("text", x = 8, y = 0.8,
           label = "4.2%", color = "firebrick")

#explore denoised counts ###################################################################################################################
denoisecount <- BASiCS_DenoisedCounts(SCEwt, Chainwt)
colSums(denoisecount)

write.csv(denoisecount,"wt_denoisecunt.csv")

megM <- as.matrix(denoisecount)

##plot PCA check batch effect
M <- t(megM)
#remove constant/zero column
Mc <- M[ , which(apply(M, 2, var) != 0)]
##PCA
pca <- prcomp(Mc, scale = TRUE)
#plot by PC1 and PC2
plot(pca$x[,3],pca$x[,2])
#the percentage of variance accounted by PC
pca.var <-  pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="screen plot",xlab = "PC",ylab = "percent variation")
##############################################################################################################################################