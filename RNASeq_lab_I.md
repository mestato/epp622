## Software installation

### Set your python version to 2.7

```{php}
module load python/2.7.11
python -V ## this should output 'Python 2.7.11'
```

### Install anaconda

```{php}
cd ~
wget https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86_64.sh
bash Anaconda2-4.1.1-Linux-x86_64.sh ## accept or confirm all prompts during installation

## re-logon to newton to activate anaconda
conda ## you should expect to see the help document for conda if you have installed anaconda successfully

conda config --add channels bioconda ## add the bioconda channel: a repository of ~1500 bioinformatic packages.
```

### Install bioinformatic tools

```{php}
conda install star -y
conda install hisat2 -y
conda install rampmap -y
conda install samtools -y
conda install skewer -y
```

## Data preparation (from *Arabidopsis thaliana*)

### Set up data directory

```{php}
cd ~ 
mkdir RNASeq_lab_I
cd RNASeq_lab_I && mkdir 0_raw_data
```

### Reference genome

```{php}
cd ~/RNASeq_lab_I/0_raw_data
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa.gz
```

### Reference transcriptome

```{php}
cd ~/RNASeq_lab_I/0_raw_data
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz
```

### Genome annotation

```{php}
cd ~/RNASeq_lab_I/0_raw_data
wget ftp://ftp.ensemblgenomes.org/pub/release-28/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.28.gtf.gz
gunzip *gz ## uncompress all .gz files in your current directory
```

### Raw reads ([Experimental info](https://github.com/mestato/epp622/blob/master/RNA_labs_data/experimental_info.csv))

```{php}
cd ~/RNASeq_lab_I/0_raw_data
cp /data/home/mchen33/EPP622_2016_fall/rnaseq_labs_data.tar.gz
tar -xvzf rnaseq_labs_data.tar.gz 
```

## RNASeq alignment/mapping: *STAR*, *HISAT2* and *RapMap*

### STAR

1. __Index reference genome__

    ```{php}
    cd ~/RNASeq_lab_I
    make alignment_STAR && cd alignment_STAR
    mkdir genomeDir
    
    star --runMode genomeGenerate     \
        --genomeDir ./genomeDir       \
        --genomeFastaFiles ../0_raw_data/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa  \
        --runThreadN 4      \
        --sjdbGTFfile ../0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf   \
        --sjdbOverhang 101
    ```
    
    * `runMode genomeGenerate`: run genome indices generation job, default is to run alignment.
    * `--genomeDir`: specify the directory for storing genome indices
    * `--genomeFastaFiles`: one or more FASTA files with genome reference sequences
    * `--runThreadN`: the number of threads to use.
    * `sjdbGTFfile`: The annotation file that STAR uses to build splice junctions database
    * `sjdbOverhang`: specifies the length of genomic sequence around the annotated junction. Usually it is set to __*Readlength - 1*__.
    
    + Command line to get the read length
        
        * Read 1: `head -2 ../0_raw_data/DRR016140_1.1percent.fastq | awk "{print length}" | tail -2`
        * Read 2: `head -2 ../0_raw_data/DRR016140_2.1percent.fastq | awk "{print length}" | tail -2`
    

2. __Align the reads__

    * We have 16 pairs of reads file. We can use a for loop to get the job done, instead of running the command 16 times.
        
        ```{php}
        cd ~/RNASeq_lab_I/alignment_STAR ## make sure you are in the right directory
        mkdir alignment_output  ## create a directory to store the alignment output files
        
        for i in `seq 25 40`
        do
             star --genomeDir ./genomeDir       \
                  --readFilesIn ../0_raw_data/DRR0161${i}_1.1percent.fastq ../0_raw_data/DRR0161${i}_2.1percent.fastq      \
                  --outFileNamePrefix ./alignment_output/DRR0161${i}_  \
                  --outSAMtype BAM SortedByCoordinate     \
                  --runThreadN 4
        done
        ```
        
        * `--genomeDir`: specifies the directory where you put your genome indices
        * `--readFilesIn`: your paired RNASeq reads files.
        * `--outFileNamePrefix`: your output file name prefix. 
        * `--outSAMtype`: your output file type. Here we want the generated bam file to be sorted by coordination.
        * `--runThreadN`: the number of threads to be used.

### HISAT2

1. __Index reference genome__

    ```{php}
    cd ~/RNASeq_lab_I
    mkdir alignment_hisat2 && cd alignment_hisat  ## create a directory for hisat2 alignment
    mkdir genomeDir  ## again, we create a directory for the genome indices
    
    hisat2-build -p 4 \
                 ../0_raw_data/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa
                 ./genomeDir/Athal_index
    ```
    
    * `-p`: specifies the number of threads to use
    * `../0_raw_data/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa`: path to the reference genome
    * `./genomeDir/Athal_index`: the base of indices files that will be generated.

2. __Align the reads__
    
    ```{php}
    cd ~/RNASeq_lab_I/alignment_hisat2 ## go into the correct directory
    mkdir alignment_output ## create a directory to store the alignment output files
    
    for i in `seq 25 40`
    do
         hisat2 ./genomeDir/Athal_index      \
                -1 ../0_raw_data/DRR0161${i}_1.1percent.fastq   \
                -2 ../0_raw_data/DRR0161${i}_2.1percent.fastq   \
                -S ./aligned_output/DRR0161${i}.sam             \
                -p 4
    done
    ```
    
    * `./genomeDir/Athal_index`: path to the directory of genome indices
    * `-1`: specifies the first paired reads file
    * `-2`: specifies the second paired reads file
    * `-S`: output to a SAM file. Default is stdout.
    * `-p`: the number of threads to use.

3. __Generate sorted bam files__

    * __hisat2__ does not generate sorted bam file automatically. We will need to use __samtools__ to convert the sam files to sorted bam fles.
    
    ```{php}
    cd ~/RNASeq_lab_I/alignment_hisat2
    mkdir sorted_bam
    
    for i in `seq 25 40`
    do
       samtools sort ./alignment_output/DRR0161${i}.sam -o ./sorted_bam/DRR0161${i}_sorted.bam
    done
    ```
    
### RapMap: *align reads to the transcriptome instead of the genome* 

1. __Index reference genome__
    
    ```{php}
    cd ~/RNASeq_lab_I
    mkdir alignment_rapmap && cd alignment_rapmap
    
    rapmap quasiindex   \
           -t ../0_raw_data/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa \
           -i ./transcriptomeDir
    ```
    
    * `quasiindex`: builds a suffix array-based index
    * `-t`: specifies the path to the reference transcriptome file
    * `i`: specifies the location where the index should be written

2. __Align the reads__
    
    ```{php}
    cd ~/RNASeq_lab_I/alignment_rapmap
    mkdir alignment_output
    
    for i in `seq 25 40`
    do
        rapmap quasimap \
               -i ./transcriptomeDir    \
               -1 ../0_raw_data/DRR0161${i}_1.1percent.fastq    \
               -2 ../0_raw_data/DRR0161${i}_2.1percent.fastq    \
               -o ./alignment_output/DRR0161${i}.sam
    done
    ```
    
    * `quansimap`: map reads using the suffix-array based method, should match the method you used for indexing.
    * `-1`: specifies the first set of reads from a paired library.
    * `-2`: specifies the second set of reads from a paired library.
    * `-o`: path to the file where the output should be written.
