## Count reads

### Count reads from STAR alignment

```{php}

for sorted_bam in $(find . -name *.bam | grep $1)
do
    htseq-count -f bam \
                -t gene \
                -s no \
                -i gene_id  \
                sorted bam \
                /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf > 
done
```

### Count reads from hisat2 alignment

### Count reads from RapMap alignment
