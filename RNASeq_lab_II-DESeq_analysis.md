## Count reads

### Write a script to automatically extract counts from all sorted bam files

Run command `cd ~/RNASeq_lab_I` and then create a script file named __`count_reads.sh`__. Put the following content into __`count_reads.sh`__. 

```{php}
for sorted_bam_path in $(find ~/RNASeq_lab_I -name *.bam | grep $1)
do
    counts_file=~/RNASeq_lab_I/counts_$1/$(echo $sorted_bam_path | grep -o "DRR0161[0-9]*")_$1_ct
    echo "The target bam file is: "$sorted_bam_path
    echo "==================================================="
    htseq-count -f bam \
                -t gene \
                -s no \
                -i gene_id  \
                $sorted_bam_path \
                /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf > $counts_file
    echo "The count data has been written into: $counts_file"
    echo "==================================================="
done
```

Change the file mode to make it an executable.

```{php}
chmod u+x count_reads.sh 
```

### Count reads from STAR alignment

### Count reads from hisat2 alignment

### Count reads from RapMap alignment
