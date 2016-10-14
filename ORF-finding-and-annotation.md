### Transcript fasta file (from trinity)

```{php}
wget 
```

### ORF finding

```{php}
TransDecoder.LongOrfs -t Trinity.fasta
```

* Outputs

```{R}
longest_orfs.pep   : all ORFs meeting the minimum length criteria, regardless of coding potential.
longest_orfs.gff3  : positions of all ORFs as found in the target transcripts
longest_orfs.cds   : the nucleotide coding sequence for all detected ORFs
```

* How many sequences in the transcript fasta file?

```{php}
grep '^>' Trinity.fasta | wc -l
```

* How many sequences have ORF?

```{php}
grep '^>' *pep | wc -l
```


