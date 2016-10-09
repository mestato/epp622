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

__Run the command line below if you are aligned to get count data from the STAR mapping results.__

```{php}
./count_reads.sh STAR
```

### Count reads from hisat2 alignment

__Run the command line below if you are aligned to get count data from the hisat mapping results.__

```{php}
./count_reads.sh hisat2
```

## Count matrix

### Preprocess data with unix command lines
[Example count data](https://github.com/mestato/epp622/blob/master/RNA_labs_data/example_count_data_hisat2.csv)
```{php}
echo gene_ID $(ls DRR* | grep -o "DRR0161[0-9]*" | tr "\n" ' ') | tr -s [:blank:] ',' > count_data.csv
paste $(ls DRR* | sort) | awk '{for(i=3;i<=NF;i+=2) $i=""}{print}' | tr -s [:blank:] ',' >> count_data.csv
```

[Experimental information](https://github.com/mestato/epp622/blob/master/RNA_labs_data/experimental_info.csv)

### Load data into R

```{R}
countData = read.csv('count_data.csv', header = TRUE, row.names = 1)
colData = read.csv("https://raw.githubusercontent.com/mestato/epp622/master/RNA_labs_data/experimental_info.csv", header = TRUE, row.names = 2)[, c("factor1", "factor2")]

## construct the data that analyzing functions from DESeq2 can recognize.
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ factor1 + factor2)
dds
```


* Pre-filtering: discard rows that have 0 for all treatments

```{R}
dim(dds)  ## before filtering
dds = dds[rowSums(counts(dds))>1, ]
dim(dds)
```

* Differential expression analysis
        + `log2 fold change (MAP): factor2 saline vs ABA` means that the estimates are log2(treated/untreated)
```{R}
dds = DESeq(dds)
res = results(dds)
res
```


