## Count reads

### Count reads from STAR alignment

```{php}

cd ~/RNASeq_lab_I
mkdir counts_$1 && cd counts_$1
for sorted_bam_path in $(find ~/RNASeq_lab_I -name *.bam | grep $1)
do
    counts_file=$(grep -o DRR0161[0-9]* $sorted_bam_path)_$1
    htseq-count -f bam \
                -t gene \
                -s no \
                -i gene_id  \
                $sorted_bam_path \
                /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf > ~/RNASeq_lab_I/counts_$1/$counts_file
done
```

### Count reads from hisat2 alignment

### Count reads from RapMap alignment
