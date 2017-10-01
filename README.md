## PAIRADISE: Paired analysis of Allelic Differential Isoform Expression

### Requirements
R version 3.2.1 or higher.

### Installation
Install the PAIRADISE .tar files then type the command ```install.packages("path_to_file/PAIRADISE", repos = NULL, type = "source")``` within R.

### Usage:
``` pairadise(my.data, numCluster, sig.level, nIter, tol, pseudocount, equal.variance) ```

The input format for the dataframe required by PAIRADISE should be as follows:

Each row of the dataframe corresponds to a different alternative splicing event. The dataframe should have 7 columns, arranged as follows: 

1. Column 1 contains the ID of the alternative splicing events. 
2. Column 2 contains counts of isoform 1 corresponding to the first group.
3. Column 3 contains counts of isoform 2 corresponding to the first group.
4. Column 4 contains counts of isoform 1 corresponding to the second group.
5. Column 5 contains counts of isoform 2 corresponding to the second group.
6. Column 6 contains the effective length of isoform 1.
7. Column 7 contains the effective length of isoform 2.

Replicates in columns 2-5 should be separated by commas, e.g. "1623,432,6" for three replicates and the replicate order should be consistent for each column to ensure pairs are matched correctly. 

Other (optional) inputs to pairadise include:

1. numCluster: Number of clusters to use for parallel computing. Default is ```numCluster = 2```.
2. sig.level: The desired level of statistical significance. Default is ```sig.level = 0.01```.
3. nIter: The maximum number of iterations of the optimization algorithm allowed. Default is ```nIter = 100```.
4. tol: Specifies the tolerance level for terminating the optimization algorithm, defined as the difference in log-likelihood ratios between iterations. Default is ```tol = 10^(-2)```.
5. pseudocount: Specifies a value for a pseudocount added to each count (e.g. values in columns 2-5 of the input dataframe) at the beginning of the analysis. Default is ```pseudocount = 0```.
6. equal.variance: Are the group variances assumed equal? Takes value TRUE or FALSE. Default is ```equal.variance = FALSE```.


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

### Example:

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
```



## Allele Specific Alternative Splicing (ASAS) Detection Pipeline 

### Dependencies & Required Installations:
Python 2.6 or higher <br>
rPGA 2.0 (https://github.com/Xinglab/rPGA) <br>
STAR (https://github.com/alexdobin/STAR)  <br>
PAIRADISE (see above sections for installation and useage)<br>

### Required Scripts (Can be downloaded from ASASDetectPipeline.tar.gz)
processGTF.SAMs.py, id2gene.py, count.py, FDR.py, PAIRADISE_fast.R, run.PAIRADISE_with_filter.sh

### Reference and Annotation Bundle
#### Option 1: Use the same genome reference and annotation used in our manuscript
https://drive.google.com/drive/folders/0B6Gm87MT7rC5WE5TMlJyQnlpdGM <br>
What's included: hg19 reference genome, Ensembl GRCh37 v75 gtf file and RADAR known RNA editing sites v2. For our specific analysis, NA12878 VCF files (1KGP phase3) are also included.
#### Option 2: Use the latest version
1. Ensembl FTP
reference genome-> http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
genome annotation -> http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
2. GATK provides a bundle of reference genome and gtf file <br>
https://software.broadinstitute.org/gatk/download/bundle
3. User can also find reference genome from UCSC genome browser and genome annotation from GENCODE 

### Toy Example
https://drive.google.com/drive/folders/0B6Gm87MT7rC5U1dDRlVNMUluVzA?usp=sharing

### Pipeline Outline
STEP I. Mapping & Preparing Allele-Specific AS Counts (using rPGA)-> 
STEP II. Statistical test (using PAIRADISE) -> 
STEP III. Result Filtering and FDR Controlling

#### STEP I. Mapping & Preparing ASAS Counts (mainly use rPGA package)
Memory Requirement: This step includes reference genome generation, STAR mapping and parsing aligned BAM files, therefore, it is recommended to allocate memory sufficiently. eg. for a sample like ENCSR000AED, 35 GB will suffice to run mapping and assigning reads, 15 GB will suffice for producing ASAS counts.
Commands: Although each rPGA command is commented in below section, users are feel free to use rPGA -h to look at the manual provided by the package. Additional useage of rPGA please refer to https://github.com/Xinglab/rPGA.

1.1 Generating the Personal Genome (hap1.fa and hap2.fa), followed by STAR indexing
```
rPGA personalize \ 'personalize' is a rPGA command to generate reference genome based known SNPs (VCF files) 
 -o /path/to/personal/genome/ \ # output directory
 -v /path/to/VCF/directory \ # eg. 'NA12878' folder in the shared 'Reference Bundle' folder
 -r /path/to/reference/genome/hg19.fa \ # can be replaced by users
 --rnaedit \ # # mask RNA editing sites
 -e /path/to/known/RNA/editing/sites/Human_AG_all_hg19_v2.txt \ # can be replaced by users
 --gz # use when input VCF files are gzipped
# MAIN OUTPUT: /path/to/personal/genome/hap1.fa,/path/to/personal/genome/hap2.fa

# Genome coordinates in VCF files should contain no 'chr'. For example use '1' instead of 'chr1'.
# A STAR indexing step is REQUIRED in order to use the next step 'rPGA mapping' command. STAR index please see STAR manual.
# Put STAR index files under the same directory with new folder HAP1 HAP2, eg. /path/to/personal/genome/HAP1 
```
1.2 Mapping to hap1 and hap2 using STAR, and assigning haplotype specific reads.
```
OUTPUT_BAM=/path/to/store/BAM/files

rPGA mapping \ 'mapping' is a command of rPGA package to run STAR to map input reads to two haplotype genomes
 -o $OUTPUT_BAM \ # output directory
 -s INPUT1.fq.gz,INPUT2.fq.gz \ # input fastq
 --hap \ # enable mapping to two personal reference genomes (hap1 and hap2)
 -g /path/to/genome/annotation/Homo_sapiens.Ensembl.GRCh37.75.gtf \ # can be replaced by users 
 --genomedir /path/to/personal/genome/HAP1/STARindex,/path/to/personal/genome/HAP2/STARindex # index generated by last step
 -N 6 \ # threads value passed to STAR
 --readlength READ_LENGTH \ # readlength value passed to STAR
 --gz \ # use when input fastq files are gzipped
# MAIN OUTPUT: $OUTPUT_BAM/HAP1/STARalign/Aligned.out.sorted.bam, $OUTPUT_BAM/HAP2/STARalign/Aligned.out.sorted.bam

rPGA assign \ # 'assign' is a command of rPGA package to assign aligned reads (generated by 'mapping') to alleles.
 -o $OUTPUT_BAM \ # output directory defined by last command ('mapping' command)
 -v /path/to/VCF/directory \ # eg. 'NA12878' folder in the shared 'Reference Bundle' folder
 -e /path/to/known/RNA/editing/sites/Human_AG_all_hg19_v2.txt \ #Can be replaced by users
 --rnaedit \ # mask RNA editing sites
 --gz # use when VCF files are gzipped
# MAIN OUTPUT: $OUTPUT_BAM/hap1.sorted.bam, $OUTPUT_BAM/hap2.sorted.bam

```
[optional] Alternaively, this assign command can also run by each chromosome and merge them later to speed up. 

1.3 Generating AS Events
```
# cd to one level higher than OUTPUT_BAM directory. 
# If you have multiple samples (multiple $OUTPUT_BAM), we suggest to finish 1.1-2 steps for all of them before run this step.
python /path/to/rMATs/processGTF.SAMs.py \
 /path/to/genome/annotation/Homo_sapiens.Ensembl.GRCh37.75.gtf \ # can be replaced by users 
 ASEvents/fromGTF \ # a standard GTF annotation output prefix
 $OUTPUT_BAM/hap1.sorted.bam,$OUTPUT_BAM/hap2.sorted.bam \ # output generated from 'assign' step
 fr-unstranded \
 temp
# MAIN OUTPUT: a set of standard GTF annotation output under ASEvents, eg fromGTF.SE.txt
```
1.4 Generating ASAS Counts file
```
# stay at the same level of directory with ASEvents folder generated last step.
cat samples.txt | while read s # samples.txt is a list of paths to BAMs generated in assign step, one path per line.  
do    
rPGA splicing \ # a funcion in rPGA to find allele specific splicing and output the junction counts
 -o ${s} \ # sample directory listed in samples.txt
 --asdir ASEvents \ # select AS events folder
 --readlength READ_LENGTH \ # read length passed to STAR
 --anchorlength ANCHOR_LENGTH # We used --readlength 100 --anchorlength 8 for ENCSR000AED
done
# MAIN OUTPUT: Sample ASAS count file, eg. ASAS.SNP.SE.JunctionReadsOnly.byPair.txt

# Merge ASAS counts of all samples to one (you need to run this even if you only have one sample in samples.txt):
ASASCounts=/path/to/ASASCounts/results
mkdir $ASASCounts
rPGA splicing \ # a funcion in rPGA to find allele specific splicing and output the junction counts
 --merge \ # merge count files for multiple samples
 --pos2id \ # change snp pos to snp id in count files
 --samples samples.txt \ # see above
 -o ASASCounts \ # output directory to store merged ASAS counts files
 -v /path/to/VCF/directory # eg. 'NA12878' folder in the 'Reference Bundle'
 # MAIN OUTPUT: Merged ASAS count file, eg. ASASCounts/ASAS.SNP.SE.JunctionReadsOnly.byPair.unfiltered.txt
 # The ourput generated this step can be directly used by PAIRADISE. See STEP II
 # We used 'unfiltered' output as we will do all filter in STEP III.
```
#### STEP II. Run PAIRADISE
Memory requirement: Less than 1 GB. <br>
Below is an example R code to run PAIRADISE. <br>
For a ready-to-use script, see 'PAIRADISE_fast.R' in ASASDetectPipeline.tar.gz. <br>
Note that we provided 'run.PAIRADISE_with_filter.sh' to perform both PAIRADISE and the post filtering steps in one single command. <br> 
```
args<-commandArgs(TRUE)
library('PAIRADISE')
my.data=read.table(args[1],header = TRUE,colClasses = c(rep("character", 5), rep("numeric", 2)))
results <- pairadise(my.data, numCluster = 8, equal.variance = F)
write.table(cbind(results$exonID,results$raw.pvalues),file=paste0(args[1],'.pairadise_out.txt'))
```

#### STEP III. Filtering for Significant ASAS Events
Memory Requirment: Less than 2 GB. <br>
Below are steps to further refining PAIRADISE raw test result for your analysis. <br>

1. Filtering low coverage events (Counts >=10)
2. Filtering small PSI values (abs(delta PSI) larger than 0.05)
3. FDR calculation and filtration (Based on FDR 10%)

We provided a simple script 'run.PAIRADISE_with_filter.sh' to perform both PAIRADISE and the post filtering steps. <br>
```
bash run.PAIRADISE_with_filter.sh script_location ASASCountDir AS_type
```
Details please see ASASDetectPipeline.tar.gz. <br>



### Contacts and bug reports:
Yi Xing: yxing@ucla.edu

Levon Demirdjian: levondem@ucla.edu

Shihao Shen: shihao@ucla.edu

If you found a bug or mistake in this project, we would like to know about it. Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been fixed.
2. Check that your input is in the correct format and you have selected the correct options.
3. Please reduce your input to the smallest possible size that still produces the bug; we will need your input data to reproduce the problem, and the smaller you can make it, the easier it will be.

