##### Download data 

```shell
#configue computing resource on Biowulf
sinteractive --cpus-per-task=32 --mem=64g --time=10:00:00 --gres=lscratch:50 --tunnel

#download contact matrix
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/2bf4cce8-6e03-422b-a190-e51d4a07d501/4DNFI1E6NJQJ.hic
#download AB compartment bigwig data
wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/a906c6b2-1563-40d4-bbde-6913c43ca8d4/4DNFIAXG3ZPL.bw


#download bigwigtobedgraph on the local laptop
rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigWigToBedGraph ./
#convert bigwig to bedgraph on the local laptop. The bedgraph can be used as a bed file
./bigWigToBedGraph /Users/lyuj2/Downloads/HAP1compartment.bw HAP1compartment.BedGraph

```

##### Call TADS

```shell
#convert hic to cool for resolution 10K,100K
hicConvertFormat -m HAP1.hic --inputFormat hic --outputFormat cool -o wt_matrix.cool --resolutions 10000 100000

#convert cool to h5
hicConvertFormat -m wt_matrix_100000.cool --inputFormat cool --outputFormat h5 -o wt_matrix_100K.h5

#check reads distribution
hicCorrectMatrix diagnostic_plot -m wt_matrix_100K.h5 -o hic_corrected.png

#find TADs
hicFindTADs -m wt_matrix_100K.h5 \
--outPrefix wt_min300k_max1M_step100k_thres0.05_delta0.01_fdr \
--minDepth 300000 \
--maxDepth 1000000 \
--step 100000 \
--thresholdComparisons 0.05 \
--delta 0.01 \
--correctForMultipleTesting fdr \
-p 32
```

##### Map genes to TADs by bedops

```shell
#activate my conda on biowulf
source /data/$USER/conda/etc/profile.d/conda.sh  
   conda activate python3.8
   
#sort cordinates for TADs(a.k.a. TADmiddle), TADboundaries, AB compartments
sort-bed TADbound_0.05.bed > TADbound_0.05.sorted.bed
sort-bed TADmiddle.bed > TADmiddle.sorted.bed
sort-bed ABcompartments.bed > ABcompartments.sort.bed

#map genes to TADs(a.k.a. TADmiddle), TADboundaries, AB compartments
bedmap --fraction-ref 1 --echo --echo-map-id-uniq --delim '\t' gene.coor.sorted.bed TADbound_0.05.sorted.bed > wt_gene_TADbound_0.05.bed
bedmap --fraction-ref 1 --echo --echo-map-id-uniq --delim '\t' gene.coor.sorted.bed TADmiddle.sorted.bed > wt_gene_TADmiddle.bed
bedmap --fraction-ref 1 --echo --echo-map-id-uniq --delim '\t' gene.coor.sorted.bed ABcompartments.sort.bed > wt_gene_ABcomp.bed
```

