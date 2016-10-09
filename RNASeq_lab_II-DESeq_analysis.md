## Get sorted bam files


## Count reads

### Write a script to automatically extract counts from all sorted bam files

Run command `cd ~/RNASeq_lab_I` and then create a script file named __`count_reads.sh`__. Put the following content into __`count_reads.sh`__. 

```{php}
#!/bin/bash

## USAGE
## ./count_reads.sh hisat
## ./count_reads.sh STAR

mkdir counts_$1 && cd counts_$1
for sorted_bam_path in $(find ./ -name *.bam | grep $1)
do
    counts_file=~/RNASeq_lab_I/counts_$1/$(echo $sorted_bam_path | grep -o "DRR0161[0-9]*")_$1_ct
    echo "The target bam file is: "$sorted_bam_path
    echo "==================================================="
    htseq-count -f bam \
                -t gene \
                -s no \
                -i gene_id  \
                $sorted_bam_path \
                /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf \
                | \
                grep -v '^__' > $counts_file
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

Change the file mode to make it an executable.

```{php}
chmod u+x count_reads.sh
```

### Count reads from STAR alignment

* __Change directory to where you have these two directory: `alignment_STAR` and `alignment_hisat2`__

        + If you use your own output, you should do `cd ~/RNASeq_lab_I`.
        + If you copy the results from my directory, you should do `cd ~/alignment_outputs`.

__Run the command line below if you want to get count data from the STAR mapping results.__

```{php}
./count_reads.sh STAR
## needs 9m56.775s
```

### Count reads from hisat2 alignment

__Run the command line below if you want to get count data from the hisat mapping results.__

```{php}
./count_reads.sh hisat2
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

* __Transfer the file `count_data.csv` to your local computer with firezilla or the `scp` command.__

    + run the `scp command` on your local computer: 
    
            = `scp newton_user@newton.login.utk.edu:/path/to/count_data.csv /path/to/directory/you/want/to/put/count_data.csv`
 
### Load data into R

Open Rstudio and change working directory to where you store the file __`count_data.csv`__.

```{R}
install.packages("DESeq2")
library("DESeq2")

countData = read.csv('count_data.csv', header = TRUE, row.names = 1)
colData = read.csv("https://raw.githubusercontent.com/mestato/epp622/master/RNA_labs_data/experimental_info.csv", header = TRUE, row.names = 2)[, c("phenotype", "stress")]

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

Results explanation:



  + `log2 fold change (MAP): factor2 saline vs ABA` means that the estimates are log2(treated/untreated)
        
```{R}
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds = DESeq(dds, test="LRT", reduced = ~phenotype)
res = results(dds)
res
```





## Interaction effect

### Method 1: create a combined group

```{R}
dds$group <- factor(paste0(dds$phenotype, dds$stress))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
```

### Method 2: likelihood ratio

```{R}
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress + phenotype:stress)
dds = DESeq(dds)
dds = DESeq(dds, test="LRT", reduced = ~phenotype+stress)
resultsNames(dds)
res = results(dds)
res$padj[order(res$padj)]
```
