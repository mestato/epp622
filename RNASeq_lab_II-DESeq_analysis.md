
### Count Data

```{php}
wget https://github.com/mestato/epp622/blob/master/RNA_labs_data/count_data.csv?raw=true
```

### Install DESeq2

```{R}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

install.packages("colorspace")
require("DESeq2")
```

## Get sorted bam files

Option 1: Use your own sorted bam files from last lab (either from hisat2 mapping or STAR mapping). 

Option 2: Copy the sorted bam files from my directory to your home directory.

```{php}
cd ~
cp -r /data/home/mchen33/RNASeq_lab_2_DESeq ./
```

## Install htseq
This requires internet access, so do this from the head node.
~~~
conda install htseq -y
~~~

## Count reads

Get an interactive session first:

```{php}
qrsh
```

### Write a script to automatically extract counts from all sorted bam files

Run command `cd ~/RNASeq_lab_2_DESeq` and then create a script file named __`count_reads.sh`__. Put the following content into __`count_reads.sh`__. 

```{php}
#!/bin/bash

## USAGE
## ./count_reads.sh hisat /path/to/the/sorted_bam_file_directory bam_file_ORDER_TYPE
## ./count_reads.sh STAR /path/to/the/sorted_bam_file_directory bam_file_ORDER_TYPE

mkdir counts_$1

for sorted_bam_path in $(find $2 -name *.bam)
do
    counts_file=$(echo $sorted_bam_path | grep -o "DRR0161[0-9]*")_$1_ct
    echo "The target bam file is: "$sorted_bam_path
    echo "==================================================="
    htseq-count -f bam \
                -t gene \
                -i gene_id  \
                -r $3 \
                $sorted_bam_path \
                ~/RNASeq_lab_2_DESeq/Arabidopsis_thaliana.TAIR10.28.gtf \
                | \
                grep -v '^__' > ./counts_$1/$counts_file
    echo "The count data has been written into: $counts_file"
    echo "==================================================="
done
```

* Why do we need this command line: `grep -v '^__'`?
```
...
...
...
ATMG01370	0
ATMG01380	0
ATMG01390	272
ATMG01400	0
ATMG01410	0
__no_feature	8798
__ambiguous	4637
__too_low_aQual	1423
__not_aligned	2794
__alignment_not_unique	52464
```

Change the file mode to make it executable.

```{php}
chmod u+x count_reads.sh 
```

### Count reads from STAR alignment



__Run the command line below if you want to get count data from the STAR mapping results.__

```{php}
./count_reads.sh STAR ~/RNASeq_lab_2_DESeq/STAR_alignment_output pos
## needs 9m56.775s
```

### Count reads from hisat2 alignment

__Run the command line below if you want to get count data from the hisat mapping results.__

```{php}
./count_reads.sh hisat2 ~/RNASeq_lab_2_DESeq/hisat2_sorted_bam name
## needs 10m4.509s
```



## Count matrix

* __Change directory to `counts_hisat2` OR `counts_STAR`__.

### Preprocess data with unix command lines
[Example count data](https://github.com/mestato/epp622/blob/master/RNA_labs_data/example_count_data_hisat2.csv)
```{php}
echo gene_ID $(ls DRR* | grep -o "DRR0161[0-9]*" | tr "\n" ' ') | tr -s [:blank:] ',' > count_data.csv
paste $(ls DRR* | sort) | awk '{for(i=3;i<=NF;i+=2) $i=""}{print}' | tr -s [:blank:] ',' >> count_data.csv
```

[Experimental information](https://github.com/mestato/epp622/blob/master/RNA_labs_data/experimental_info.csv)

* __Transfer the file `count_data.csv` to your local computer with filezilla or the `scp` command.__

    + run the `scp command` on your local computer: 
    
    ```
    scp newton_user@newton.login.utk.edu:/path/to/count_data.csv /path/to/directory/you/want/to/put/count_data.csv
    ```
 
### Load data into R

Open Rstudio and change working directory to where you store the file __`count_data.csv`__.

```{R}
install.packages("DESeq2")
library("DESeq2")

countData = read.csv('count_data.csv', header = TRUE, row.names = 1)
colData = read.csv("https://raw.githubusercontent.com/mestato/epp622/master/RNA_labs_data/experimental_info.csv", header = TRUE, row.names = 2)[, c("phenotype", "stress")]
colnames(countData) = paste0(colData$phenotype, '_', colData$stress)
rownames(colData) = paste0(colData$phenotype, '_', colData$stress)

## construct the data that analyzing functions from DESeq2 can recognize.
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds
```


* Pre-filtering: discard rows that have 0 for all treatments

```{R}
dim(dds)  ## before filtering
dds = dds[rowSums(counts(dds))>1, ]
dim(dds)
```

### Differential expression analysis

```{R}
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds = dds[rowSums(counts(dds))>1, ]                             
dds = DESeq(dds)
res = results(dds)
res
```

Results

```{R}
log2 fold change (MAP): stress saline vs ABA 
Wald test p-value: stress saline vs ABA 
DataFrame with 19453 rows and 6 columns
             baseMean log2FoldChange     lfcSE        stat       pvalue       padj
            <numeric>      <numeric> <numeric>   <numeric>    <numeric>  <numeric>
AT1G01010   0.9246572     0.04106919 1.1664635  0.03520829   0.97191365         NA
AT1G01020   0.7345512    -0.50228525 1.1997038 -0.41867437   0.67545413         NA
AT1G01030   0.4918562     0.53399011 1.2159565  0.43915231   0.66055118         NA
AT1G01040   4.1881473    -0.88294298 0.7686287 -1.14872497   0.25066941  0.5074783
AT1G01050  12.9203449     0.95553494 0.4842900  1.97306342   0.04848834  0.1769671
...               ...            ...       ...         ...          ...        ...
ATMG01350   0.1237857     -0.6264714 0.9607553  -0.6520613 0.5143615962         NA
ATMG01360   0.4759696      1.1442281 1.2182994   0.9392011 0.3476275109         NA
ATMG01370   0.1996883     -0.4199930 1.0363781  -0.4052508 0.6852931724         NA
ATMG01380   0.1976165      0.3782464 1.0249796   0.3690282 0.7121066974         NA
ATMG01390 235.3293034      0.8171322 0.2399232   3.4058071 0.0006596878 0.00617547
```


__Results explanation:__

We choose gene ATMG01390 as an example

+ __baseMean__: the average of normalized counts across all samples. This represents the intercept of your GLM.
+ __log2Foldchange__: (stress saline vs ABA): log2(treated/untreated). Here it is log2(saline/ABA)
+ __lfcSE__: standard error of the log2FoldChange estimate
+ __stat__: statistic for the hypothesis test. For Wald test, it is log2FoldChange/lfcSE. 
+ __pvalue__: the corresponding p-value from Wald test (or likelihood ratio test)
+ __padj__: adjusted p value due to multiple comparisons.
+ __reasons for NA values__:
   
    * This gene has 0 count for all samples
    * This row has an extreme count outlier
    * Filtered by automatic independent filtering for having low mean count.
   
   
Get a summary of the results:

```{R}
summary(res)
```

By default, the `results()` function display comparison of the last level of the last variable over the first level of this variable. 

We can extract results of comparisons between other levels:

```{R}
results(dds, contrast=c("stress", "mock", "saline"))
```

We can also extract results of comparisons from the *phenotype* factor

```{R}
results(dds, contrast=c("phenotype", "ros1-3", "wildtype"))
```

__MA-plot of the results__

```{R}
plotMA(res)
```

__Plot counts__

We select the gene which has the smallest padj value to plot its count at different levels

```{R}
mostSigGene = rownames(dds)[which.min(res$padj)]
mostSigGene
par(mfcol=c(1,2))
plotCounts(dds, gene=mostSigGene, intgroup="stress")
plotCounts(dds, gene=mostSigGene, intgroup="phenotype")
par(mfcol=c(1,1))
```

### Exporting results into CSV files

```{R}
orderedRes = res[order(res$padj), ]
sigOrderedRes = subset(orderedRes, padj < 0.05)
write.csv(as.data.frame(sigOrderedRes), "STRESS_saline_vs_ABA.csv")
```

### Likelihood Ratio Test (LRT) method

The LRT method compares the full model with the reduced model to see if the *removed variable (factor)* has significant effect on the fitted model.

```{R}
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds = DESeq(dds, test="LRT", reduced = ~phenotype)
resLRT = results(dds)
resLRT
```

### Wald test vs Likelihood Ratio test

How many genes are significant in LRT test?

```{R}
orderedResLRT = resLRT[order(resLRT$padj), ]
sigOrderedResLRT = subset(orderedResLRT, padj<0.05)
dim(sigOrderedResLRT)
```

How many genes are significant in Wald test?

```{R}
orderedRes = res[order(res$padj), ]
sigOrderedRes = subset(orderedRes, padj < 0.05)
dim(sigOrderedRes)
```

__Why is the number from LRT much larger than the number from Wald test?__

Let's a gene that are significant in LRT but not in Wald test.

```{R}
LRTsigGenes = rownames(sigOrderedResLRT)
Res_genesSigInLRT = res[LRTsigGenes, ]
notSigWaldGene = rownames(Res_genesSigInLRT)[which.max(Res_genesSigInLRT$padj)]
notSigWaldGene
```

padj comparison for the selected gene
```{R}
resLRT[notSigWaldGene, ]
res[notSigWaldGene, ]
```

Let's plot the counts for this gene

```{R}
plotCounts(dds, gene=notSigWaldGene, intgroup="stress")
```

What do you see from the plot?

Let's extract the results of comparisons between __dehydration__ and "ABA".

```{R}
res2 = results(dds, contrast=c("stress", "dehydration", "ABA"))
res2[notSigWaldGene, ]
```


