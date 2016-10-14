### Create a working directory

```{php}
mkdir 18_annotation
cd 18_annotation/
```

### Transcript fasta file (output from trinity assembly)

```{php}
wget https://raw.githubusercontent.com/mestato/epp622/master/RNA_labs_data/Trinity.fasta
```

### ORF finding

* Installation of  transdecoder
    + It seems newton's transdecoder does not work. It needs some perl module to be installed. So we are going to install it with conda.

        ```{R}
        conda install transdecoder -y
        TransDecoder.LongOrfs
        ```

* Run
    ```{php}
    TransDecoder.LongOrfs -t Trinity.fasta
    ```

* Outputs
    
    ```{php}
    cd Trinity.fasta.transdecoder_dir/
    ls
    ```
    
    ```{R}
    longest_orfs.pep   : all ORFs meeting the minimum length criteria, regardless of coding potential.
    longest_orfs.gff3  : positions of all ORFs as found in the target transcripts
    longest_orfs.cds   : the nucleotide coding sequence for all detected ORFs
    ```

* How many sequences have ORF?

    ```{php}
    grep '^>' *pep | wc -l
    ```

### Functional annotation

* __Pfam search__: search the peptides for protein domains using Pfam.

    + Get software and Pfam database ready
    ```{php}
    PATH=$PATH:/lustre/projects/rnaseq_ws/apps/hmmer-3.1b2-linux-intel-x86_64/binaries
    
    mkdir Pfam_search  ## create a directory for the Pfam database
    cd Pfam_search
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm ## prepare an HMM database for faster hmmscan searches
    ```

    + Run
    ```{php}
    hmmscan -o my_hmmscan.out   \
            --tblout my_hmmscan.SeqHits.tblr    \
            --domtblout my_hmmscan.DomainHits.tblr  \ 
            -E 1e-5 \
            ./Pfam-A.hmm ../longest_orfs.pep
    ```
    
* __blastP search__: search a protein database

    + Get software and protein database ready
    
    ```{php}
    module load blast
    module switch blast
    
    mkdir uniprot_blastp_search
    cd uniprot_blastp_search
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot
    ```
    
    + Run
    
    ```{php}
    blastp  -query longest_orfs.pep  \
            -db uniprot_sprot.fasta  
            -max_target_seqs 1 
            -outfmt 6 
            -evalue 1e-5 
            -num_threads 10 > blastp.outfmt6
    ```


### Improved ORF finding

__*We can include homology searches as ORF retention criteria*__


