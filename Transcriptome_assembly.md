### Set Up

We will need slightly more memory to run a de novo assembly, so lets request an interactive session with 7Gb of RAM

```{php}
 qrsh -l mem=7G -q medium*
```


Navigate to our working directory. Lets create a folder for this analysis

```{php}
cd ~
mkdir trainscriptome_trinity_assembly
cd trainscriptome_trinity_assembly
```

Get raw RNAseq data

```{php}
cp -r /data/home/mchen33/EPP622_2016_fall/0_raw_reads .
```

We will be using the assembler Trinity.


### Run

Install __*trinity*__

```{php}
conda install trinity -y
Trinity --help
```

There are a number of required parameters

```{php}
--seqType <string>  :type of reads: ( fa, or fq )
--max_memory <string>  :suggested max memory to use by Trinity where limiting can be enabled.
```

If paired reads:

```{php}
--left  <string>    :left reads, one or more file names (separated by commas, no spaces)
--right <string>    :right reads, one or more file names (separated by commas, no spaces)
```

Or, if unpaired reads:

```{php}
--single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
```

So we will need to get all the forward reads in one file, and all the reverse reads in another. So lets concatenate the raw reads files

```{php}
 cat 0_raw_reads/*1.1percent.fastq > allR1.fastq
 cat 0_raw_reads/*2.1percent.fastq> allR2.fastq
```

Other parameters of use?

```{php}
--SS_lib_type <string>  :Strand-specific RNA-Seq read orientation.
--min_contig_length <int>  :minimum assembled contig length to report
--CPU <int>  :number of CPUs to use, default: 2
--jaccard_clip high gene density, potential UTR overlap
```

Trinity needs to be able to find the software programs bowtie and samtools.

```{php}
 PATH=$PATH:/lustre/projects/rnaseq_ws/apps/bowtie-1.1.1
 PATH=$PATH:/lustre/projects/rnaseq_ws/apps/samtools-1.1/
```

The default java version is 1.8. We need to load java 1.7 for trinity.

```{php}
module load java/jre7u60
java -version ## check java version
```

Now we have figured out all our parameters, so lets run the assembly software. We will give it 6Gb of RAM instead of 7Gb so that it does not use too much and kill the interactive session.

```{php}
Trinity \
 --seqType fq \
 --max_memory 6G \
 --left allR1.fastq \
 --right allR2.fastq \
 --SS_lib_type FR \
 --CPU 1 \
 --min_contig_length 60 \
 --jaccard_clip
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
