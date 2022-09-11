###### Experimental Design

```
1, Modeling count data with Negative Binomial distribution
2, To minimize bench effect, equally assign confounders to each group.
3, How many replicates are enough? the more the better, typically 3. I'd like to have biological 4 replicates since it significantly increases the number of detected DE genes.
4, sequencing depth. typically 30 M per samples
```

###### QC for raw reads

```shell
fastqc
cutadapt
Qualimap2

A tool called Qualimap explores the features of aligned reads in the context of the genomic region they map to, hence providing an overall view of the data quality (as an HTML file). Various quality metrics assessed by Qualimap include:

DNA or rRNA contamination
5’-3’ biases
Coverage biases
```

###### Reads mapping

```shell
#generate index for 50bp mapping
STAR --runThreadN 16 --runMode genomeGenerate\
 --genomeDir /data/lyuj2/ref/index/STAR/GRCH38_v35/\
 --genomeFastaFiles /data/lyuj2/ref/index/STAR/GRCH38_v35/GRCh38.p13.genomeE.fa\
 --sjdbGTFfile /data/lyuj2/ref/index/STAR/GRCH38_v35/gencode.v35.annotationE.gtf\
 --sjdbOverhang 49
 
#mapping
 STAR --runThreadN 16\
     --genomeDir /data/lyuj2/ref/index/STAR/GRCH38_v35/\
     --readFilesIn UMI_R_F.fastq\
     --outFileNamePrefix ./UMI_R\
     --outFilterMultimapNmax 1\
     --outSAMtype BAM SortedByCoordinate
     
     
#or use salmon, the standarded way nowadays.
```

###### Counting

```shell
#featurecounts, count reads mapping to genes/exons
GTF=/Users/lyuj2/Desktop/RNAData/Reference/Homo_sapiens.GRCh38.101.gtf
 featureCounts -a $GTF -o RNAtp_A_exon_counts.txt RNAtp_A.bam
 
We are going to use the following options:

-T 4 # specify 4 cores
-p species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
-s specifies strand-specific read counting. 0 for unstranded reads, 1 for stranded reads and 2 for reversely stranded reads. This depends on the library used in the sequencing protocol.

and the following are the values for the required parameters:

-a ~/unix_lesson/rnaseq/reference_data/chr1-hg19_genes.gtf # required option for specifying path to GTF

-o ~/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt # required option for specifying path to, and name of the text output (count matrix)

~/unix_lesson/rnaseq/results/STAR/bams/*bam # the list of all the bam files we want to collect count information for

# or use salmon, the standarded way nowadays
```

###### QC based on counts (DEseq2)

```R
1, Create a DESeqDataSet object, dds

2, normalize data by DESeq2’s median of ratios, which considers sequencing depth and RNA composition. 
The gene count comparisons between samples and for DE analysis; NOT for within sample comparisons

3, get normalized and rlog counts for PCA and Hierarchical Clustering
############################################################################################################
#Create DESEq2 object
the DESEq2 object is a pre-specified list, containing a count matrix and a metadata. The data stored in these pre-specified slots can be accessed by using specific package-defined functions.We will also need to specify a design formula. The design formula specifies the column(s) in the metadata table and how they should be used in the analysis.

#take data from the salmon result
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

#generate size factors and store it 
dds <- estimateSizeFactors(dds)

#check the size. usually these size factors are around 1. If you see large variations between samples, it is important to take note as it might indicate the presence of extreme outliers.
sizeFactors(dds)

#normalize and retrieve the normalized counts matrix from dds
normalized_counts <- counts(dds, normalized=TRUE)
#save this normalized data matrix to file for later use
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#######################################################################################################
visulize the data for QC
#log normalized count matrix using rlog and vst
NOTE1: “Many common statistical methods for exploratory analysis of multidimensional data, especially methods for clustering and ordination (e. g., principal-component analysis and the like), work best for (at least approximately) homoskedastic data; this means that the variance of an observable quantity (i.e., here, the expression strength of a gene) does not depend on the mean. In RNA-Seq data, however, variance grows with the mean. For example, if one performs PCA directly on a matrix of normalized read counts, the result typically depends only on the few most strongly expressed genes because they show the largest absolute differences between samples. A simple and often used strategy to avoid this is to take the logarithm of the normalized count values plus a small pseudocount; however, now the genes with low counts tend to dominate the results because, due to the strong Poisson noise inherent to small count values, they show the strongest relative differences between samples.

As a solution, DESeq2 offers the regularized-logarithm transformation, or rlog for short. For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation. For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples. Using an empirical Bayesian prior in the form of a ridge penality, this is done such that the rlog-transformed data are approximately homoskedastic.” - From the “Beginner’s guide to using the DESeq2 package” by Love, Anders and Huber, 2014 (the DESeq2 vignette is the updated version of this doc).

NOTE2: The DESeq2 vignette suggests large datasets (100s of samples) to use the variance-stabilizing transformation (vst) instead of rlog for transformation of the counts, since the rlog function might take too long to run and the vst() function is faster with similar properties to rlog.

# Transform counts
rld <- rlog(dds, blind=TRUE)
The blind=TRUE argument is to make sure that the rlog() function does not take our sample groups into account - i.e. does the transformation in an unbiased manner.The rlog() function returns a DESeqTransform object, another type of DESeq-specific object. The reason you don’t just get a matrix of transformed values is because all of the parameters (i.e. size factors) that went into computing the rlog transform are stored in that object. We use this object to plot the PCA and heirarchical clustering figures for quality assessment.
NOTE: The rlog() funtion can be a bit slow when you have e.g. > 20 samples. In these situations the vst() function is much faster and performs a similar transformation appropriate for use with plotPCA().
### Plot PCA 
plotPCA(rld, intgroup="sampletype")
### Extract the rlog matrix from the object
rld_mat <- assay(rld)
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
### Load pheatmap package
library(pheatmap)
### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = meta)
```

###### DE analysis

```R
if QC is fine, The final step in the differential expression analysis workflow is fitting the raw counts to the NB model and performing the statistical test for differentially expressed genes. A design formula tells the statistical software the known sources of variation to control for, as well as, the factor of interest to test for during differential expression testing. For example, if you know that sex is a significant source of variation in your data, then sex should be included in your model. The design formula should have all of the factors in your metadata that account for major sources of variation in your data. The last factor entered in the formula should be the condition of interest.If you want to examine the expression differences between treatments, and you know that major sources of variation include sex and age, then your design formula would be:

design = ~ sex + age + treatment

The tilde (~) should always precede your factors and tells DESeq2 to model the counts using the following formula. Note the factors included in the design formula need to match the column names in the metadata.

## Run analysis
dds <- DESeq(dds)

#QC for DE analysis
Plot dispersion estimates to examine to ensure your data is a good fit for the DESeq2 model. You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing mean expression levels.

plotDispEsts(dds)


```

###### Functional analysis (GO)

```R
# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)

require Ensembl and Entrez IDs.
#res_ids, your dataset with Ensembl or Entrez IDs
## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(res_ids$gene)

## Extract significant results. If very few genes are in any of these lists (< 50, roughly) it may not be possible to get any significant GO terms.
sigOE <- dplyr::filter(res_ids, padj < 0.05)

sigOE_genes <- as.character(sigOE$gene)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")

#################visulization########################################
## Dotplot 
The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 GO terms by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
  
dotplot(ego, showCategory=50)


####the enrichment GO plot
which shows the relationship between the top 50 most significantly enriched GO terms (padj.), by grouping similar terms together. Before creating the plot, we will need to obtain the similarity between terms using the pairwise_termsim() function (instructions for emapplot). In the enrichment plot, the color represents the p-values relative to the other displayed terms (brighter red is more significant), and the size of the terms represents the number of genes that are significant from our list.
## Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50)

##the category netplot 
the category netplot shows the relationships between the genes associated with the top five most significant GO terms and the fold changes of the significant genes associated with these terms (color). The size of the GO terms reflects the pvalues of the terms, with the more significant terms being larger. This plot is particularly useful for hypothesis generation in identifying genes that may be important to several of the most affected processes.

Note - You may need to install the ggnewscale package using install.packages("ggnewscale") for the cnetplot() function to work.

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange

names(OE_foldchanges) <- sigOE$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)
         
## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)

cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)
```

###### Functional analysis (GSEA)

```R
The hypothesis of FCS methods is that although large changes in individual genes can have significant effects on pathways (and will be detected via ORA methods), weaker but coordinated changes in sets of functionally related genes (i.e., pathways) can also have significant effects. Thus, rather than setting an arbitrary threshold to identify ‘significant genes’, all genes are considered in the analysis. The gene-level statistics from the dataset are aggregated to generate a single pathway-level statistic and statistical significance of each pathway is reported. This type of analysis can be particularly helpful if the differential expression analysis only outputs a small list of significant DE genes.

1.to perform the analysis, we will need to acquire the Entrez IDs.

2.we need to order the fold changes in decreasing order. To do this we’ll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse’s arrange(), which requires a data frame.

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

######################################Performing GSEA###########################################################
First, we will set the seed so that we all obtain the same result:

set.seed(123456)
NOTE: The permutations are performed using random reordering, so every time we run the function we will get slightly different results. If we would like to use the same permutations every time we run a function, then we use the set.seed(123456) function prior to running. The input to set.seed() can be any number, but if you would want the same results, then you would need to use the same number as the lesson.

To perform the GSEA using KEGG gene sets with clusterProfiler, we can use the gseKEGG() function:

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
              organism = "hsa", # supported organisms listed below
              nPerm = 1000, # default number permutations
              minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
              pvalueCutoff = 0.05, # padj cutoff value
              verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
NOTE: The nPerm argument was left at its default value of 1000. This parameter specifies how many times this randomization (for perumtations) is done. The more randomizations that are performed, the more precise the FDR q-value estimation will be individual terms/pathways. Therefore, if you are finding few or no terms enriched you might want to try increasing this number.
  
Explore the GSEA plot of enrichment of one of the pathways in the ranked list:

## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03008')
                                                                                                                            
 Use the Pathview R package to integrate the KEGG pathway data from clusterProfiler into pathway images:

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
              pathway.id = "hsa03008",
              species = "hsa",
              limit = list(gene = 2, # value gives the max/min limit for foldchanges
              cpd = 1))
NOTE: If the below error message occurs: Error in detach("package:dplyr", unload = T) : invalid 'name' argument, that means the dplyr package is not currently loaded. Ignore the message and continue to run pathview command.

NOTE: Printing out Pathview images for all significant pathways can be easily performed as follows:

## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
   pathview(gene.data = foldchanges, 
            pathway.id = gseaKEGG_results$ID[x], 
            species = "hsa",
            limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), 
           get_kegg_plots)
There are other gene sets available for GSEA analysis in clusterProfiler (Disease Ontology, Reactome pathways, etc.). In addition, it is possible to supply your own gene set GMT file, such as a GMT for MSigDB using special clusterProfiler functions as shown below:
  
BiocManager::install("GSEABase")
library(GSEABase)

# Load in GMT file of gene sets (we downloaded from the Broad Institute for MSigDB)

c2 <- read.gmt("/data/c2.cp.v6.0.entrez.gmt.txt")

msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)

msig_df <- data.frame(msig)
 
                                                                                                                            
                                                                                                                                                                                                                                                      Resources for functional analysis
g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi
DAVID - http://david.abcc.ncifcrf.gov/tools.jsp
clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
GeneMANIA - http://www.genemania.org/
GenePattern - http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
AmiGO - http://amigo.geneontology.org/amigo
ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/
WGCNA - https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
GSEA - http://software.broadinstitute.org/gsea/index.jsp
SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
GAGE/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html
```

###### GeneID conversion

```R

Annotation databases :
OrgDb, TxDb, Go.db, EnsDb, and BioMart annotations
e.g.
org.Hs.eg.db
There are a plethora of organism-specific orgDb packages, such as org.Hs.eg.db for human and org.Mm.eg.db for mouse, and a list of organism databases can be found here. These databases are best for converting gene IDs or obtaining GO information for current genome builds, but not for older genome builds. These packages provide the current builds corresponding to the release date of the package, and update every 6 months. If a package is not available for your organism of interest, you can create your own using AnnotationHub.


Interface tool: 
AnnotationDbi: queries the OrgDb, TxDb, Go.db, EnsDb, and BioMart annotations.
AnnotationDbi is an R package that provides an interface for connecting and querying various annotation databases using SQLite data storage. The AnnotationDbi packages can query the OrgDb, TxDb, EnsDb, Go.db, and BioMart annotations. 

  
 # Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Check object metadata
org.Hs.eg.db

We can easily extract information from this database using AnnotationDbi with the methods: columns, keys, keytypes, and select. For example, we will use our org.Hs.eg.db database to acquire information, but know that the same methods work for the TxDb, Go.db, EnsDb, and BioMart annotations.

# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                     keys = res_tableOE_tb$gene,  # data to use for retrieval
                                     columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retreive for given data
                                     keytype = "ENSEMBL") # type of data given in 'keys' argument

 Let’s take a peek to see if we actually returned annotations for each individual Ensembl gene ID that went in to the query:

length(which(is.na(annotations_orgDb$SYMBOL)))

Let’s get rid of those NA entries:

# Determine the indices for the non-NA genes
non_na_idx <- which(is.na(annotations_orgDb$SYMBOL) == FALSE)

# Return only the genes with annotations using indices
annotations_orgDb <- annotations_orgDb[non_na_idx, ]

de-duplicate our data before using it.

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]




##############################generate the Ensembl annotations#####################################

To generate the Ensembl annotations, the EnsDb database can also be easily queried using AnnotationDbi. You will need to decide the release of Ensembl you would like to query. We know that our data is for GRCh38, and the most current EnsDb release for GRCh38 in Bioconductor is release 86, so we can install this database. All Ensembl releases are listed here. NOTE: this is not the most current release of GRCh38 in the Ensembl database, but it’s as current as we can obtain through AnnotationDbi.

Since we are using AnnotationDbi to query the database, we can use the same functions that we used previously:

# Load the library
library(EnsDb.Hsapiens.v86)

# Check object metadata
EnsDb.Hsapiens.v86

# Explore the fields that can be used as keys
keytypes(EnsDb.Hsapiens.v86)
Now we can return all gene IDs for our gene list:

# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys = res_tableOE_tb$gene,
                                           columns = c("SYMBOL", "ENTREZID","GENEBIOTYPE"),
                                           keytype = "GENEID")
We can check for NA entries, and find that there are none:

length(which(is.na(annotations_edb$SYMBOL) == FALSE))
Then we can again deduplicate, to remove the gene symbols which appear more than once:

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_edb <- annotations_edb[non_duplicates_idx, ]
NOTE: In this case we used the same build but a slightly older release, and we found little discrepancy. If your analysis was conducted using an older genome build (i.e hg19), but used a newer build for annotation some genes may be found to be not annotated (NA). Some of the genes have changed names in between versions (due to updates and patches), so may not be present in the newer version of the database.
```



###### Workflow (Summary)

```R
Summary of differential expression analysis workflow
We have detailed the various steps in a differential expression analysis workflow, providing theory with example code. To provide a more succinct reference for the code needed to run a DGE analysis, we have summarized the steps in an analysis below:

Obtaining gene-level counts from Salmon using tximport

 # Run tximport
 txi <- tximport(files, 
         type="salmon", 
         tx2gene=t2g, 
         countsFromAbundance = "lengthScaledTPM")
	
 # "files" is a vector wherein each element is the path to the salmon quant.sf file, and each element is named with the name of the sample.
 # "t2g" is a 2 column data frame which contains transcript IDs mapped to geneIDs (in that order)
Creating the dds object:

 # Check that the row names of the metadata equal the column names of the **raw counts** data
 all(colnames(txi$counts) == rownames(metadata))
	
 # Create DESeq2Dataset object
 dds <- DESeqDataSetFromTximport(txi, 
                 colData = metadata, 
                 design = ~ condition)
Exploratory data analysis (PCA & hierarchical clustering) - identifying outliers and sources of variation in the data:

 # Transform counts for data visualization
 rld <- rlog(dds, 
         blind=TRUE)
	
 # Plot PCA 
 plotPCA(rld, 
     intgroup="condition")
	
 # Extract the rlog matrix from the object and compute pairwise correlation values
 rld_mat <- assay(rld)
 rld_cor <- cor(rld_mat)
	
 # Plot heatmap
 pheatmap(rld_cor, 
      annotation = metadata)
Run DESeq2:

     # **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation using "dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ covaraite + condition)"

 # Run DESeq2 differential expression analysis
 dds <- DESeq(dds)

     # **Optional step** - Output normalized counts to save as a file to access outside RStudio using "normalized_counts <- counts(dds, normalized=TRUE)"
Check the fit of the dispersion estimates:

 # Plot dispersion estimates
 plotDispEsts(dds)
Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:

 # Specify contrast for comparison of interest
 contrast <- c("condition", "level_to_compare", "base_level")
	
 # Output results of Wald test for contrast
 res <- results(dds, 
            contrast = contrast, 
            alpha = 0.05)
	
 # Shrink the log2 fold changes to be more accurate
 res <- lfcShrink(dds, 
          coef = "sampletype_group1_vs_group2", 
          type = "apeglm")	 
      # The coef will be dependent on what your contras was. and should be identical to what is stored in resultsNames()
Output significant results:

 # Set thresholds
 padj.cutoff < - 0.05
	
 # Turn the results object into a tibble for use with tidyverse functions
 res_tbl <- res %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>% 
               as_tibble()
	
 # Subset the significant results
 sig_res <- filter(res_tbl, 
           padj < padj.cutoff)
Visualize results: volcano plots, heatmaps, normalized counts plots of top genes, etc.

Perform analysis to extract functional significance of results: GO or KEGG enrichment, GSEA, etc.

Make sure to output the versions of all tools used in the DE analysis:

sessionInfo()
For better reproducibility, it can help to create RMarkdown reports, which save all code, results, and visualizations as nicely formatted html reports. We have available an example html report for perusal. To create these reports we have additional materials available.
```

