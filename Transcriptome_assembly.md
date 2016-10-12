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

We will be using the assembler Trinity.

Run
Lets find out more information about the software and its parameters

 ls /lustre/projects/rnaseq_ws/apps/trinityrnaseq_r20131110/
 /lustre/projects/rnaseq_ws/apps/trinityrnaseq_r20131110/Trinity.pl --help
There are a number of required parameters

--seqType <string>  :type of reads: ( fa, or fq )
--JM <string>  :(Jellyfish Memory) number of GB of system memory to use for
If paired reads:

--left <string>  :left reads, one or more (separated by space)
--right <string>  :right reads, one or more (separated by space)
Or, if unpaired reads:

--single <string>  :single reads, one or more (note, if single file contains pairs, can use flag: --run_as_paired )
So we will need to get all the forward reads in one file, and all the reverse reads in another. So lets concatenate the raw reads files

 cat /lustre/projects/rnaseq_ws/raw_data/Sp*R1.fq > allR1.fastq
 cat /lustre/projects/rnaseq_ws/raw_data/Sp*R2.fq > allR2.fastq
Other parameters of use?

--SS_lib_type <string>  :Strand-specific RNA-Seq read orientation.
--min_contig_length <int>  :minimum assembled contig length to report
--CPU <int>  :number of CPUs to use, default: 2
--jaccard_clip high gene density, potential UTR overlap

Trinity needs to be able to find the software programs bowtie and samtools.

 PATH=$PATH:/lustre/projects/rnaseq_ws/apps/bowtie-1.1.1
 PATH=$PATH:/lustre/projects/rnaseq_ws/apps/samtools-1.1/

Now we have figured out all our parameters, so lets run the assembly software. We will give it 6Gb of RAM instead of 7Gb so that it does not use too much and kill the interactive session.

 /lustre/projects/rnaseq_ws/apps/trinityrnaseq_r20131110/Trinity.pl \
 --seqType fq \
 --JM 6G \
 --left allR1.fastq \
 --right allR2.fastq \
 --SS_lib_type FR \
 --CPU 1 \
 --min_contig_length 60 \
 --jaccard_clip
