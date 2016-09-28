## Count reads

### Count reads from STAR alignment

```{php}
htseq-count -f bam \
            -t gene \
            -s no \
            -i gene_id  \
            /data/home/mchen33/RNASeq_lab_I/alignment_STAR/alignment_output/DRR016125_Aligned.sortedByCoord.out.bam \
            /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf
```

### Count reads from hisat2 alignment

### Count reads from RapMap alignment
