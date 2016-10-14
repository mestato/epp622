### Transcript fasta file (from trinity)

```{php}
wget 
```

### ORF finding

```{php}
TransDecoder.LongOrfs -t Trinity.fasta
```

* Outputs

```{php}

```

* How many sequences in the transcript fasta file?

```{php}
grep '^>' Trinity.fasta | wc -l
```

* How many sequences have ORF?

```{php}
grep '^>' *pep | wc -l
```


