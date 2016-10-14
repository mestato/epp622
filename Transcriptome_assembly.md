### Set Up

We will need slightly more memory to run a de novo assembly, so lets request an interactive session with 7Gb of RAM

```{php}
 qrsh -l mem=7G -q medium*
```


Navigate to our working directory. Lets create a folder for this analysis

```{php}
cd ~
mkdir 16_trinity_assembly
cd 16_trinity_assembly
```

Get raw RNAseq data

```{php}
cp -r /data/home/mchen33/EPP622_2016_fall/0_raw_reads/DRR016125* .
cp -r /data/home/mchen33/EPP622_2016_fall/0_raw_reads/DRR016126* .
```

We will be using the assembler Trinity.


### Run

__*Trinity*__ is already installed. We just need to load it.

```{php}
module load trinity/2.2.0
```

We also need samtools and bowtie2 to run trinity assembly.

```{php}
module switch samtools/1.3.1 ## samtools is already loaded, so we need to use switch command to switch to a different version
module load bowtie2/2.2.8
```

* There are a number of required parameters
  
  ```{php}
  --seqType <string>  :type of reads: ( fa, or fq )
  --max_memory <string>  :suggested max memory to use by Trinity where limiting can be enabled.
  ```

* If paired reads:
  
  ```{php}
  --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
  --right <string>    :right reads, one or more file names (separated by commas, no spaces)
  ```

* Or, if unpaired reads:
  
  ```{php}
  --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
  ```


* Other parameters of use?
  
  ```{php}
  --SS_lib_type <string>  :Strand-specific RNA-Seq read orientation.
  --min_contig_length <int>  :minimum assembled contig length to report
  --CPU <int>  :number of CPUs to use, default: 2
  --jaccard_clip high gene density, potential UTR overlap
  ```



The default java version is 1.8. We need to load java 1.7 for trinity.

```{php}
module load java/jre7u60
java -version ## check java version
```

Now we have figured out all our parameters, so lets run the assembly software. We will give it 6Gb of RAM instead of 7Gb so that it does not use too much and kill the interactive session. Create a script file named `trinity_assembly.qsh` which has the following code in it.

```{php}
#$ -N trinity_assembly
#$ -cwd
#$ -S /bin/bash
#$ -l mem=7G
#$ -q medium*

module switch samtools/1.3.1 ## samtools is already loaded, so we need to use switch command to switch to a different version
module load bowtie2/2.2.8
module load java/jre7u60 
module load trinity/2.2.0

Trinity \
 --seqType fq \
 --max_memory 6G \
 --left DRR016125_1.1percent.fastq,DRR016126_1.1percent.fastq \
 --right DRR016125_2.1percent.fastq,DRR016126_2.1percent.fastq \
 --SS_lib_type FR \
 --CPU 1 \
 --min_contig_length 60 \
 --jaccard_clip
```

Submit the job:

```{php}
qsub trinity_assembly.qsh
```

### Results

All the output files are placed into a directory named trinity_out_dir (you can change this in the parameters if you want)

The assembled transcripts can be found in

```{php}
  ./trinity_out_dir/Trinity.fasta
```

A quick one line command to check the number of seqeunces in the file

```{php}
 grep -c '^>' Trinity.fasta
```

Most of the output files can be ignored. Lets keep looking at Trinity.fasta. Meg wrote a perl script that can report statistics about the file.

```{php}
 /lustre/projects/rnaseq_ws/apps/fasta_file_stats.pl Trinity.fasta
```

### Assessing quality of assembly with Transrate

```{php}
PATH=$PATH:/lustre/projects/rnaseq_ws/apps/transrate-1.0.1-linux-x86_64/
```

```{php}
transrate --assembly trinity_out_dir/Trinity.fasta \
          --left DRR016125_1.1percent.fastq,DRR016126_1.1percent.fastq 
          --right DRR016125_2.1percent.fastq,DRR016126_2.1percent.fastq > transrate_output
```
