# PAIRADISE: Detecting allele-specific alternative splicing from population-scale RNA-seq data

## About

We introduce PAIRADISE (PAIred Replicate analysis of Allelic DIfferential Splicing Events), a method for detecting allele-specific alternative splicing (ASAS) from RNA-seq data. PAIRADISE uses a statistical model that aggregates ASAS signals across multiple individuals in a population. It formulates ASAS detection as a statistical problem for identifying differential alternative splicing from RNA-seq data with paired replicates. The PAIRADISE statistical model is applicable to many forms of allele-specific isoform variation (e.g. RNA editing), and can be used as a generic statistical model for RNA-seq studies involving paired replicates.

## Installation

Install the dependencies:
- Python >= (2.7.3)
- R >= (3.2.1)
- rPGA >= (2.0.0)
- STAR >= (2.4.2a)
- pysam >= (0.8.4)
- pybedtools >= (0.7.4)

Download PAIRADISE from github:
	git clone github.com/Xinglab/PAIRADISE

Install the PAIRADISE statistical model in R from local repository:
	install.packages('PAIRADISE_1.0.tar.gz', repos = NULL)

## Input files
Test examples of the input files are available at TestExample.tar.gz
- RNA-seq Data (Example: chr20.test.1.fastq.gz,chr20.test.2.fastq.gz)
- VCF Genotype Files (Example: 20.vcf.gz)

## Annotation Files

- hg19 reference genome: http://hgdownload.cse.ucsc.edu/downloads.html#human
- Ensemble GRCh37 v75 GTF file: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
- RADAR known RNA editing sites v2: http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt

## Test PAIRADISE:

Download and unzip the PAIRADISE test dataset and annotation files: https://sourceforge.net/projects/pairadise/files/TestExample.tar.gz/download

    tar xvfz TestExample.tar.gz
    python PAIRADISE.py -o TEST_chr20 -v TestExample/20.vcf.gz -r TestExample/chr20.hg19.fa -e Human_AG_all_hg19_v2.txt -gtf Homo_sapiens.Ensembl.GRCh37.75.gtf -s TestExample/chr20.test.1.fastq.gz,TestExample/chr20.test.2.fastq.gz -gz

## Usage

    python PAIRADISE.py -h
    
    usage: usage: PAIRADISE.py [options] arg1 arg2
	
    optional arguments:
		  -h, -help            show this help message and exit
		  -version             Version.
		  -gtf GTF             An annotation of genes and transcripts in GTF format.
		  -r ref               FASTQ file of reference genome.
		  -e rnaedit           List of RNA editing sites.
		  -s fastq             FASTQ files of RNA-seq data, multiple individuals separated by comma.
		  -v VCF               Folders of VCF files, multiple individuals separated by comma.
		  -o OD                Output folder of post step.
		  -N Number            Max number of read mismatches, default is 3
		  -M Number            Max number of multiple alignments, default is 20. 
		  -readLength LENGTH   The read length.
		  -anchorLength LENGTH Anchor length of junction read, default is 8.		 
		  -statoff             Turn statistical analysis off.
		  -gz                  Use gzipped FASTQ and VCF files.

## Output

All output files are in outputFolder:

- AS_Event.MATS.JC.txt evaluates splicing with only reads that span splicing junctions
	- IC_ALLELE_1: inclusion counts for ALLELE_1, replicates are separated by comma
	- SC_ALLELE_1: skipping counts for ALLELE_1, replicates are separated by comma
	- IC_ALLELE_2: inclusion counts for ALLELE_2, replicates are separated by comma
	- SC_ALLELE_2: skipping counts for ALLELE_2, replicates are separated by comma

- AS_Event.MATS.JCEC.txt evaluates splicing with reads that span splicing junctions and reads on target (striped regions on home page figure)
	- IC_ALLELE_1: inclusion counts for ALLELE_1, replicates are separated by comma
	- SC_ALLELE_1: skipping counts for ALLELE_1, replicates are separated by comma
	- IC_ALLELE_2: inclusion counts for ALLELE_2, replicates are separated by comma
	- SC_ALLELE_2: skipping counts for ALLELE_2, replicates are separated by comma

- Important columns contained in output files above.
	- IncFormLen: length of inclusion form, used for normalization
	- SkipFormLen: length of skipping form, used for normalization
	- P-Value: the statistical significance of allele specific alternative splicing
	- FDR: the false discovery rate of allele specific alternative splicing
	- IncLevel1: inclusion level for ALLELE_1 replicates (comma separated) calculated from normalized counts
	- IncLevel2: inclusion level for ALLELE_2 replicates (comma separated) calculated from normalized counts
	- IncLevelDifference: average(IncLevel1) - average(IncLevel2)

- fromGTF.AS_Event.txt: all possible alternative splicing (AS) events derived from GTF and RNA.

## Use the Stand-Alone PAIRADISE Statistical Model in R:
PAIRADISE can be used as a stand-alone statistical model in R by downloading the and installing the file 'PAIRADISE_RPackage_1.0.tar.gz'. The following command will install the package within R:

```install.packages('PAIRADISE_RPackage_1.0.tar.gz', repos = NULL, type = 'source') ```

The primary function of the PAIRADISE R package, ```pairadise```, is called as follows:

``` pairadise(my.data) ```

Each row of the dataframe ```my.data``` corresponds to a different alternative splicing event. ```my.data``` should have 7 columns, arranged as follows: 

1. Column 1 contains the ID of the alternative splicing events. 
2. Column 2 contains counts of isoform 1 corresponding to the first group.
3. Column 3 contains counts of isoform 2 corresponding to the first group.
4. Column 4 contains counts of isoform 1 corresponding to the second group.
5. Column 5 contains counts of isoform 2 corresponding to the second group.
6. Column 6 contains the effective length of isoform 1.
7. Column 7 contains the effective length of isoform 2.

Replicates in columns 2-5 should be separated by commas, e.g. "1623,432,6" for three replicates and the replicate order should be consistent for each column to ensure pairs are matched correctly. The PAIRADISE R package includes a sample dataset satisfying the formatting requirements mentioned above; this dataset can be accessed using the command

``data("sample_dataset")``

Other (optional) inputs to ```pairadise``` include:

1. numCluster: Number of clusters to use for parallel computing. Default is ```numCluster = 2```.
2. sig.level: The desired level of statistical significance. Default is ```sig.level = 0.01```.
3. nIter: The maximum number of iterations of the optimization algorithm allowed. Default is ```nIter = 100```.
4. tol: Specifies the tolerance level for terminating the optimization algorithm, defined as the difference in log-likelihood ratios between iterations. Default is ```tol = 10^(-2)```.
5. pseudocount: Specifies a value for a pseudocount added to each count (e.g. values in columns 2-5 of the input dataframe) at the beginning of the analysis. Default is ```pseudocount = 0```.
6. seed: Seed value for random number generation.
7. equal.variance: Are the group variances assumed equal? Takes value TRUE or FALSE. Default is ```equal.variance = FALSE```.

Output:
The function ```pairadise``` returns a list containing the following entries:


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

## Example of Using the Stand-Alone PAIRADISE Statistical Model in R:

Run the following example to test if PAIRADISE is working properly:
```
## Load the PAIRADISE package
library(PAIRADISE)

## Load and save sample dataset
data("sample_dataset")
my.data <- sample_dataset

## Look at the raw data
my.data

## Run PAIRADISE and store results
results <- pairadise(my.data, equal.variance = FALSE)

