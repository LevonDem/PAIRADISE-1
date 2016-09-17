## PAIRADISE: Paired analysis of Allelic Differential Isoform Expression

### Requirements
1. Install R version 3.2.1 or higher.
2. Install the PAIRADISE package from within R, e.g. run ```install.packages("PAIRADISE")```.

### Installation
The source code can directly be called from R.

### Usage:
``` pairadise(my.data) ```

The input format for the dataframe "my.data" is described below in the Examples section.

### Examples:

Run the following example to test if PAIRADISE is working properly:
```
set.seed(12345)
nEvents <- 3  # number of alternative splicing events

## Organize data into the data frame my.data following the proper formatting:
eventID <- paste("Event", as.character(seq(1:nEvents)))
my.data <- data.frame(matrix(nrow = nEvents, ncol = 7))
my.data[,1] <- eventID
my.data[,2] <- c("12,3,5", "2,9,10,6,5,4", "15,17000,20,100")
my.data[,3] <- c("0,1,2", "0,0,4,0,3,2", "2,12,1,1")
my.data[,4] <- c("2,4,5", "12,13,7,7,7,8", "1,6,7,10")
my.data[,5] <- c("0,1,3", "0,0,0,4,3,1", "274,NA,320,5650")
my.data[,6] <- c(3,3,3)
my.data[,7] <- c(1,1,1)

## Store results
results <- pairadise(my.data)
```

The input format for the dataframe required by PAIRADISE should be in the following format:

Each row of the dataframe corresponds to a different alternative splicing event. The dataframe should have 7 columns, arranged as follows: 

1. Column 1 contains the ID of the alternative splicing events. 
2. Column 2 contains counts of isoform 1 corresponding to the first group.
3. Column 3 contains counts of isoform 2 corresponding to the first group.
4. Column 4 contains counts of isoform 1 corresponding to the second group.
5. Column 5 contains counts of isoform 2 corresponding to the second group.
6. Column 6 contains the effective length of isoform 1.
7. Column 7 contains the effective length of isoform 2.

Replicates in columns 2-5 should be separated by commas, e.g. 1623,432,6 for three replicates and the replicate order should be consistent for each column to ensure pairs are matched correctly. 

Other (optional) inputs to pairadise include:

1. numCluster: Number of clusters to use for parallel computing. Default is numCluster = 2.
2. sig.level: The desired level of statistical significance. Default is sig.level = 0.01.
3. nIter: The maximum number of iterations of the optimization algorithm allowed. Default is nIter = 100.
4. tol: Specifies the tolerance level for terminating the optimization algorithm, defined as the difference in log-likelihood ratios between iterations. Default is tol = 10^(-2).
5. pseudocount: Specifies a value for a pseudocount added to each count (e.g. values in columns 2-5 of the input dataframe) at the beginning of the analysis. Default is pseudocount = 0.


Output:
The function "pairadise" returns a list containing the following entries:

1. sig.results.Bonferroni: Matrix containing the significant exons (after Bonferroni correction at sig.level), their p-values, and test-statistics.
2. sig.results.FDR: Matrix containing the significant exons (after FDR correction using BH at sig.level), their p-values, and test-statistics.
3. testStats: Vector of test statistics for paired analysis.
4. raw.pvalues: Vector of pvalues for each exon/event.
5. param.unconstrained: List of parameter estimates for unconstrained model.
6. param.constrained: List of parameter estimates for constrained model.
7. latent.u: List of parameter estimates of latent variables for unconstrained model.
8. latent.c: List of parameter estimates of latent variables for constrained model.
9. nReplicates: Vector containing the number of valid replicates for each AS event in my.data.
10. totalIter: Vector containing total number of iterations required for optimization algorithm.
11. exonID: Character vector containing event IDs.
12. nExon: Total number of valid AS events in my.data.
13. I1: List containing all exon inclusion counts for group 1.
14. S1: List containing all exon skipping counts for group 1.
15. I2: List containing all exon inclusion counts for group 2.
16. S2: List containing all exon skipping counts for group 2.


### Contacts and bug reports:
Yi Xing: yxing@ucla.edu

Levon Demirdjian: levondem@ucla.edu

If you found a bug or mistake in this project, we would like to know about it. Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been fixed.
2. Check that your input is in the correct format and you have selected the correct options.
3. Please reduce your input to the smallest possible size that still produces the bug; we will need your input data to reproduce the problem, and the smaller you can make it, the easier it will be.

