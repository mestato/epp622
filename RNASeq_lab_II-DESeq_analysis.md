## Count reads

### Write a script to automatically extract counts from all sorted bam files

Run command `cd ~/RNASeq_lab_I` and then create a script file named __`count_reads.sh`__. Put the following content into __`count_reads.sh`__. 

```{php}
#!/bin/bash

## USAGE
## ./count_reads.sh hisat
## ./count_reads.sh STAR
## ./count_reads.sh rapmap

cd ~/RNASeq_lab_I
mkdir counts_$1 && cd counts_$1
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
                /data/home/mchen33/RNASeq_lab_I/0_raw_data/Arabidopsis_thaliana.TAIR10.28.gtf \
                | \
                grep -v '^__' > $counts_file
    echo "The count data has been written into: $counts_file"
    echo "==================================================="
done
```

* Why do we need this command line: `grep -v '^__'`?
```
...
...
...
ATMG01370	0
ATMG01380	0
ATMG01390	272
ATMG01400	0
ATMG01410	0
__no_feature	8798
__ambiguous	4637
__too_low_aQual	1423
__not_aligned	2794
__alignment_not_unique	52464
```

Change the file mode to make it an executable.

```{php}
chmod u+x count_reads.sh 
```

### Count reads from STAR alignment

__Run the command line below if you are aligned to get count data from the STAR mapping results.__

```{php}
./count_reads.sh STAR
```

### Count reads from hisat2 alignment

__Run the command line below if you are aligned to get count data from the hisat mapping results.__

```{php}
./count_reads.sh hisat2
```

## Count matrix

### Preprocess data with unix command lines
[Example count data](https://github.com/mestato/epp622/blob/master/RNA_labs_data/example_count_data_hisat2.csv)
```{php}
echo gene_ID $(ls DRR* | grep -o "DRR0161[0-9]*" | tr "\n" ' ') | tr -s [:blank:] ',' > count_data.csv
paste $(ls DRR* | sort) | awk '{for(i=3;i<=NF;i+=2) $i=""}{print}' | tr -s [:blank:] ',' >> count_data.csv
```

[Experimental information](https://github.com/mestato/epp622/blob/master/RNA_labs_data/experimental_info.csv)

### Load data into R

```{R}
countData = read.csv('count_data.csv', header = TRUE, row.names = 1)
colData = read.csv("https://raw.githubusercontent.com/mestato/epp622/master/RNA_labs_data/experimental_info.csv", header = TRUE, row.names = 2)[, c("factor1", "factor2")]

## construct the data that analyzing functions from DESeq2 can recognize.
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ factor1 + factor2)
dds
```


* Pre-filtering: discard rows that have 0 for all treatments

```{R}
dim(dds)  ## before filtering
dds = dds[rowSums(counts(dds))>1, ]
dim(dds)
```

* Differential expression analysis
        + `log2 fold change (MAP): factor2 saline vs ABA` means that the estimates are log2(treated/untreated)
```{R}
dds = DESeq(dds)
res = results(dds)
res
```


```{R}
function (object, modelMatrix = NULL, modelFormula, alpha_hat, 
    lambda, renameCols = TRUE, betaTol = 1e-08, maxit = 100, 
    useOptim = TRUE, useQR = TRUE, forceOptim = FALSE, warnNonposVar = TRUE) 
{
    if (missing(modelFormula)) {
        modelFormula <- design(object)
    }
    if (is.null(modelMatrix)) {
        modelAsFormula <- TRUE
        modelMatrix <- model.matrix(modelFormula, data = colData(object))
    }
    else {
        modelAsFormula <- FALSE
    }
    stopifnot(all(colSums(abs(modelMatrix)) > 0))
    modelMatrixNames <- colnames(modelMatrix)
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
    modelMatrixNames <- make.names(modelMatrixNames)
    if (renameCols) {
        convertNames <- renameModelMatrixColumns(colData(object), 
            modelFormula)
        convertNames <- convertNames[convertNames$from %in% modelMatrixNames, 
            , drop = FALSE]
        modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
    }
    colnames(modelMatrix) <- modelMatrixNames
    normalizationFactors <- if (!is.null(normalizationFactors(object))) {
        normalizationFactors(object)
    }
    else {
        matrix(rep(sizeFactors(object), each = nrow(object)), 
            ncol = ncol(object))
    }
    if (missing(alpha_hat)) {
        alpha_hat <- dispersions(object)
    }
    if (length(alpha_hat) != nrow(object)) {
        stop("alpha_hat needs to be the same length as nrows(object)")
    }
    if (missing(lambda)) {
        lambda <- rep(1e-06, ncol(modelMatrix))
    }
    justIntercept <- if (modelAsFormula) {
        modelFormula == formula(~1)
    }
    else {
        ncol(modelMatrix) == 1 & all(modelMatrix == 1)
    }
    if (justIntercept & all(lambda <= 1e-06)) {
        alpha <- alpha_hat
        betaConv <- rep(TRUE, nrow(object))
        betaIter <- rep(1, nrow(object))
        betaMatrix <- matrix(log2(mcols(object)$baseMean), ncol = 1)
        mu <- normalizationFactors * as.numeric(2^betaMatrix)
        logLike <- rowSums(dnbinom(counts(object), mu = mu, size = 1/alpha, 
            log = TRUE))
        deviance <- -2 * logLike
        modelMatrix <- model.matrix(~1, colData(object))
        colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
        w <- (mu^-1 + alpha)^-1
        xtwx <- rowSums(w)
        sigma <- xtwx^-1
        betaSE <- matrix(log2(exp(1)) * sqrt(sigma), ncol = 1)
        hat_diagonals <- w * xtwx^-1
        res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix, 
            betaSE = betaSE, mu = mu, betaIter = betaIter, deviance = deviance, 
            modelMatrix = modelMatrix, nterms = 1, hat_diagonals = hat_diagonals)
        return(res)
    }
    qrx <- qr(modelMatrix)
    if (qrx$rank == ncol(modelMatrix)) {
        Q <- qr.Q(qrx)
        R <- qr.R(qrx)
        y <- t(log(counts(object, normalized = TRUE) + 0.1))
        beta_mat <- t(solve(R, t(Q) %*% y))
    }
    else {
        if ("Intercept" %in% modelMatrixNames) {
            beta_mat <- matrix(0, ncol = ncol(modelMatrix), nrow = nrow(object))
            logBaseMean <- log(rowMeans(counts(object, normalized = TRUE)))
            beta_mat[, which(modelMatrixNames == "Intercept")] <- logBaseMean
        }
        else {
            beta_mat <- matrix(1, ncol = ncol(modelMatrix), nrow = nrow(object))
        }
    }
    lambdaLogScale <- lambda/log(2)^2
    betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix, 
        nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat, 
        beta_matSEXP = beta_mat, lambdaSEXP = lambdaLogScale, 
        tolSEXP = betaTol, maxitSEXP = maxit, useQRSEXP = useQR)
    mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
    dispersionVector <- rep(dispersions(object), times = ncol(object))
    logLike <- nbinomLogLike(counts(object), mu, dispersions(object))
    rowStable <- apply(betaRes$beta_mat, 1, function(row) sum(is.na(row))) == 
        0
    rowVarPositive <- apply(betaRes$beta_var_mat, 1, function(row) sum(row <= 
        0)) == 0
    betaConv <- betaRes$iter < maxit
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat
    colnames(betaMatrix) <- modelMatrixNames
    colnames(modelMatrix) <- modelMatrixNames
    betaSE <- log2(exp(1)) * sqrt(pmax(betaRes$beta_var_mat, 
        0))
    colnames(betaSE) <- paste0("SE_", modelMatrixNames)
    rowsForOptim <- if (useOptim) {
        which(!betaConv | !rowStable | !rowVarPositive)
    }
    else {
        which(!rowStable | !rowVarPositive)
    }
    if (forceOptim) {
        rowsForOptim <- seq_along(betaConv)
    }
    if (length(rowsForOptim) > 0) {
        resOptim <- fitNbinomGLMsOptim(object, modelMatrix, lambda, 
            rowsForOptim, rowStable, normalizationFactors, alpha_hat, 
            betaMatrix, betaSE, betaConv, beta_mat, mu, logLike)
        betaMatrix <- resOptim$betaMatrix
        betaSE <- resOptim$betaSE
        betaConv <- resOptim$betaConv
        mu <- resOptim$mu
        logLike <- resOptim$logLike
    }
    stopifnot(!any(is.na(betaSE)))
    nNonposVar <- sum(rowSums(betaSE == 0) > 0)
    if (warnNonposVar & nNonposVar > 0) 
        warning(nNonposVar, "rows had non-positive estimates of variance for coefficients")
    list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix, 
        betaSE = betaSE, mu = mu, betaIter = betaRes$iter, deviance = betaRes$deviance, 
        modelMatrix = modelMatrix, nterms = ncol(modelMatrix), 
        hat_diagonals = betaRes$hat_diagonals)
}
<environment: namespace:DESeq2>
```
