### Working directories

Create some directories for our analysis.

```{php}
cd ~
mkdir lesson_16
cd lesson_16
mkdir raw_data
mkdir trinity_assembly
```

Get raw RNAseq data

```{php}
cp -r /data/home/mchen33/EPP622_2016_fall/0_raw_reads/DRR016125* ./raw_data
```

Create soft links within the analysis directory

```{php}
cd trinity_assembly/
ln -s ../raw_data/* .
```

We will be using the assembler Trinity.


### Run

__*Trinity*__ is already installed. We just need to load it.

```{php}
module load trinity/2.2.0
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

Now we have figured out all our parameters, so lets run the assembly software. We will give it 6Gb of RAM instead of 7Gb so that it does not use too much and kill the interactive session. Create a script file named `trinity_assembly.qsh` which has the following code in it.

```{php}
#$ -N trinity_assembly
#$ -cwd
#$ -S /bin/bash
#$ -l mem=7G
#$ -q medium*


module load bowtie2/2.2.8
module load java/jre7u60 
module load trinity/2.2.0
module switch samtools/1.3.1 ## samtools is already loaded, so we need to use switch command to switch to a different version

Trinity \
 --seqType fq \
 --max_memory 6G \
 --left DRR016125_1.1percent.fastq \
 --right DRR016125_2.1percent.fastq \
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

The trinity assembly will take hours to finish even with very small datasets. We are going to use the results I generated a few days ago.

```{php}
cp -r /data/home/mchen33/RNASeq_trinity_assembly/trinity_out_dir_copy ./  ## copy trinity results
```

A quick one line command to check the number of seqeunces in the file

```{php}
grep -c '^>' ./trinity_out_dir_copy/Trinity.fasta
```

Most of the output files can be ignored. Lets keep looking at Trinity.fasta. Meg wrote a perl script that can report statistics about the file.

```{php}
 /lustre/projects/rnaseq_ws/apps/fasta_file_stats.pl Trinity.fasta
```

### Assessing quality of assembly with Transrate

* Get software ready

  ```{php}
  module load transrate 
  ```

  
* Run
  
  ```{php}
  transrate --assembly trinity_out_dir_copy/Trinity.fasta \
            --left DRR016125_1.1percent.fastq  \
            --right DRR016125_2.1percent.fastq > transrate_output
  ```
  
* Results

  ```{php}
  cd transrate_results
  cat assemblies.csv
  ```
  
* Modify results to make it more readable
  + option 1: copy and paste it into excel or google sheet
  + option 2: write a script
  ```{php}
  head -1 assemblies.csv | awk 'BEGIN {FS=",";}{for(i=1;i<=NF;i+=1){print $i}}' > field_name.txt
  tail -1 assemblies.csv | awk 'BEGIN {FS=",";}{for(i=1;i<=NF;i+=1){print $i}}' > field_value.txt
  paste field_name.txt field_value.txt | column -s $'\t' -t
  ```

