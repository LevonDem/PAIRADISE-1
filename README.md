# PAIRADISE: Paired Analysis of Allelic Differential Splicing Events

## About

PAIRADISE is a statistical model and computational approach to detect allele specific alternative splicing events from RNA-Seq and genotyping data. The PAIRADISE statistical model calculates the P-value and false discovery rate of the significance level that two alleles differ in alternative splicing.

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

## Annotation Files (Yang, please add the official download link to these)

- hg19 reference genome (link)
- Ensemble GRCh37 v75 GTF file (link)
- RADAR known RNA editing sites v2 (link)

## Test PAIRADISE:

Download and unzip the PAIRADISE test dataset. (Yang, please upload the test data files onto our MIMG server. Talk with Zhijie if you are not sure about how to upload.)

    tar xvfz GM12878.tar.gz
    python PAIRADISE.py -o OUTPUT/ -s1 test.chr10.1.fq.gz -s2 test.chr10.2.fq.gz -v GM12878 -r hg19.fa -gtf Homo_sapiens.Ensembl.GRCh37.75.gtf -e Human_AG_all_hg19_v2.txt -readlength 100 -gz

## Usage

    python PAIRADISE.py -h
    
    usage: usage: PAIRADISE.py [options] arg1 arg2
	
    optional arguments:
		  -h, -help            show this help message and exit
		  -version             Version.
		  -gtf GTF             An annotation of genes and transcripts in GTF format.
		  -r ref               FASTQ file of reference genome.
		  -e rnaedit           List of RNA editing sites.
		  -s1 S1               FASTQ files of the first allele, multiple individuals separated by comma.
		  -s2 S2               FASTQ files of the second allele, multiple individuals separated by comma.
		  -v VCF               Folders of VCF files, multiple individuals separated by comma.
		  -o OD                Output folder of post step.
		  -readLength LENGTH   The length of each read.
		  --statoff             Turn statistical analysis off.

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

