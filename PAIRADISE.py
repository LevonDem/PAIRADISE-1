#
## this program does deep seq data analysis
#

### import necessary libraries
import re,os,sys,logging,time,datetime;

### checking out the python version
if sys.version_info < (2,6):
  print ("Python Version error: must use phthon 2.6 or greater (not supporting python 3 yet)");
  sys.exit(-1);
#else:
#  print "Passed version test";
#  sys.exit(0);

import commands;

#
###PAIRADISE version
#
PAIRADISE_ver="1.0.0"
#
############################ parameter variables
### seven required values

#inputs
vcf=''; ## directory for VCF files, can be gzipped
fastq=''; ## rna-seq data of fastq

#annotations
gtf=''; ## gtf file, Homo_sapiens.Ensembl.GRCh37.75.gtf
ref=''; ## reference genome (FASTA), hg19.fa
edit=''; ## Human_AG_all_hg19_v2.txt

#output folder
outDir=''; ## output directory

### with default values
T=8; #number of threads STAR uses, default is 8
M=20; #max number of multiple alignments, default is 20
N=6; #max number of read pair mismatches, default is 6
readlength=75;
anchorLength=8;
gz=1; #whether the vcf and fastq are gzipped
statoff=1; #whether to run the statistical model
h=0; #whether to display tutorial
version=0; #whether to display version number

### checking out the argument names
validArgList=['-s','-v','-o','-gtf','-r','-e','-N','-M','-readLength','-anchorLength','-statoff','-gz','-h','-version'];
for argIndex in range(1,len(sys.argv)): ## going through the all parameters 
  if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
    print ('Not valid argument: %s' % sys.argv[argIndex]);
    print ('Please provide valid arguments.');
    sys.exit();

for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
  if(sys.argv[paramIndex] == '-s'):  ## rna-seq sample
    paramIndex += 1;  ## increase index
    fastq = sys.argv[paramIndex];    
  elif (sys.argv[paramIndex] == '-v'):  ## VCF folder
    paramIndex += 1;  ## increase index
    vcf = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-o'):  ## output folder
    paramIndex += 1;  ## increase index
    outDir = sys.argv[paramIndex];	
  elif (sys.argv[paramIndex] == '-gtf'):  ## gtf file
    paramIndex += 1;  ## increase index
    gtf = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-r'):  ## reference genome
    paramIndex += 1;  ## increase index
    ref = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-e'):  ## list of RNA editing 
    paramIndex += 1;  ## increase index
    edit = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-T'):  ## number of read pair mismatches
    paramIndex += 1;  ## increase index
    T = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-N'):  ## number of read pair mismatches
    paramIndex += 1;  ## increase index
    N = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-M'):  ## number of multiple alignments
    paramIndex += 1;  ## increase index
    M = int(sys.argv[paramIndex]);	
  elif (sys.argv[paramIndex] == '-readlength'):  ## read length
    paramIndex += 1;  ## increase index
    readLength = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-anchorLength'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    anchorLength = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-gz'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    gz = 1;	
  elif (sys.argv[paramIndex] == '-statoff'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    statoff = 1;
  elif (sys.argv[paramIndex] == '-h'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    h = 1;			
  elif (sys.argv[paramIndex] == '-help'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    h = 1;			
  elif (sys.argv[paramIndex] == '-version'):  ## anchor length for tophat
    paramIndex += 1;  ## increase index
    version = 1;			

if (h == 1):
  print ('usage: PAIRADISE.py [options] arg');
  print ('optional arguments:');
  print ('  		  -h, -help            show this help message and exit');
  print ('  		  -version             Version.');
  print ('  		  -gtf GTF             An annotation of genes and transcripts in GTF format.');
  print ('  		  -r ref               FASTQ file of reference genome.');
  print ('  		  -e rnaedit           List of RNA editing sites.');
  print ('  		  -s fastq             FASTQ files of RNA-seq data, multiple individuals separated by comma.');
  print ('  		  -v VCF               Folders of VCF files, multiple individuals separated by comma.');
  print ('  		  -o OD                Output folder of post step.');
  print ('  		  -N Number            Max number of read mismatches, default is 3.');
  print ('  		  -M Number            Max number of multiple alignments, default is 20.');
  print ('  		  -readLength LENGTH   The read length.');
  print ('  		  -anchorLength LENGTH Anchor length of junction read, default is 8.');
  print ('  		  -statoff             Turn statistical analysis off.');
  print ('  		  -gz                  Use gzipped FASTQ and VCF files.');
  sys.exit();

if (version == 1):
  print ('PAIRADISE 1.0.0');
  sys.exit();
 
### checking out the required arguments
if (vcf=='' or  fastq=='' or  gtf=='' or ref=='' or edit=='' or outDir==''): ### at least one required param is missing
  print ('Not enough arguments!!');
  sys.exit();
  
########## functions here... ############

def personalize(): ## get the personalized genome
  vcf.list = re.findall('[^,]+',vcf); ##assuming each RNA-seq data has its VCF file
  for vcfIndex in range(1,len(vcf.list)):
    cmd = 'mkdir '+outDir+'/'+vcfIndex;
    commands.getstatusoutput(cmd);
    cmd = 'rPGA personalize -o '+outDir+'/'+vcfIndex+'/ -v '+vcf.list[vcfIndex]+' -r '+ref+' -e '+edit+' --rnaedit';
    if gz==1:
      cmd += ' --gz ';
    commands.getstatusoutput(cmd);
  return;

def indexing(): ##Yang, please fill in the STAR indexing step
  return;
  
def mapping():
  vcf.list = re.findall('[^,]+',vcf);
  fastq.list = re.findall('[^,]+',fastq);
  for vcfIndex in range(1,len(vcf.list)):
    cmd = 'rPGA mapping -o '+outDir+'/'+vcfIndex+'/ -s '+fastq.list[vcfIndex]+' --genomedir '+outDir+'/'+vcfIndex+'/HAP1/STARindex,'+outDir+'/'+vcfIndex+'/HAP2/STARindex, -g '+gtf+' -N '+N+' -readLength '+readlength+' --hap';
    if gz==1:
      cmd += ' --gz ';
    commands.getstatusoutput(cmd);
  return;
  
def assign():
  vcf.list = re.findall('[^,]+',vcf);
  for vcfIndex in range(1,len(vcf.list)):
    cmd = 'rPGA assign -o '+outDir+'/'+vcfIndex+'/ -v '+vcf.list[vcfIndex]+' -e '+edit+' --rnaedit';
    if gz==1:
      cmd += ' --gz ';
    commands.getstatusoutput(cmd);
  return;

## bam files from all samples are used in this step
def ASASEvent():
  vcf.list = re.findall('[^,]+',vcf);
  bam.all = '';
  for vcfIndex in range(1,len(vcf.list)):
    bam.all += outDir+'/'+vcfIndex+'/hap1.sorted.bam,'+outDir+'/'+vcfIndex+'/hap2.sorted.bam,';
  bam.all = bam.all[:-1];
  cmd = 'python '+binPath+'/processGTF.SAMs.py '+gtf+' ASEvents/fromGTF '+bam.all+' fr-unstranded '+outDir+'/ASEvent/';
  commands.getstatusoutput(cmd);
  return;

## Please double check the input and output of this step. It doesn't seem right from the github.
def ASASCount():
  vcf.list = re.findall('[^,]+',vcf);
  for vcfIndex in range(1,len(vcf.list)):
    cmd = 'rPGA splicing -o '+outDir+'/'+vcfIndex+' --asdir '+outDir+'/ASEvent/ -readLength '+readlength+' -anchorLength '+anchorlength;
    commands.getstatusoutput(cmd);
  return;
 
def ASASMerge():
  vcf.list = re.findall('[^,]+',vcf);
  vcf.list.output = '';
  for vcfIndex in range(1,len(vcf.list)): ## Output the sample list for bam files
	vcf.list.output += outDir+'/'+vcfIndex+'\n';
  output = open(outDir+'/samples.txt','w');
  output.write(vcf.list.output); output.close();
  cmd = 'rPGA splicing --merge --pos2id --samples '+outDir+'/samples.txt -o '+outDir+'/ASASCounts/ -v '+vcf.list[vcfIndex];
  commands.getstatusoutput(cmd);
  return;

def PAIRADISE():
  AS.list = os.listdir(outDir+'/ASASCounts/'); ##get the list of all ASAS events
  for ASIndex in range(1, len(AS.list)):
    cmd = 'Rscript '+binPath+'/PAIRADISE_fast.R '+outDir+'/ASASCounts/'+AS.list[ASIndex];
    commands.getstatusoutput(cmd);
  return;

## Main procedure
## Debug and reporting code can be added between steps
personalize();
indexing();
mapping();
assign();
ASASEvent();
ASASCount();
ASASMerge();
PAIRADISE();
