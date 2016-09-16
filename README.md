## PAIRADISE: Paired analysis of Allelic Differential Isoform Expression

### Requirements
1. Install R version 3.2.1 or higher.
2. Install the PAIRADISE package from within R, e.g. run ```install.packages("PAIRADISE")```.

### Installation
The source code can directly be called from R.

### Usage:
Examples of code usage
``` pairadise(my.data) ```
The input format for the dataframe "my.data" is described in the Example section

### Examples:
The input format for the dataframe required by PAIRADISE should be in the following format:
Each row of the dataframe corresponds to a different alternative splicing event. The dataframe should have 7 columns, arranged as follows: Column 1 contains the ID of the exons/events. Column 2 contains counts of isoform 1 corresponding to the first group. Column 3 contains counts of isoform 2 corresponding to the first group. Column 4 contains counts of isoform 1 corresponding to the second group. Column 5 contains counts of isoform 2 corresponding to the second group. Replicates in columns 2-5 should be separated by commas, e.g. 1623,432,6 for three replicates and the replicate order should be consistent for each column to ensure pairs are matched correctly. Column 6 contains the effective length of isoform 1. Column 7 contains the effective length of isoform 2. 
Example:
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
9. nReplicates: Vector containing the number of valid replicates for each exon in my.data.
10. totalIter: Vector containing total number of iterations required for optimization algorithm.
11. exonID: Character vector containing exonIDs .
12. nExon: Total number of valid exons in my.data.
13. I1: List containing all exon inclusion counts for group 1.
14. S1: List containing all exon skipping counts for group 1.
15. I2: List containing all exon inclusion counts for group 2.
16. S2: List containing all exon skipping counts for group 2.

1. This is a list.
2. This is a list
