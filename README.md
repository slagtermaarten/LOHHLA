This is a fork of the original LOHHLA repository by Maarten Slagter, PhD student at the 
Netherlands Cancer Institute. Edits are intended to facilitate usage and interpretation of 
the program and make it more robust to various problems, while keeping the core 
functionality intact.

Currently, the major additions include:

* GRCh38 support. Do check your reference genome for presence of any additional 
  alternative HLA sequences that are not hard-coded in the source code and please let me 
  know if you find any. Perhaps these should be made user-definable.
* Removed the silent assumption of paired-end sequencing data, single-end sequencing is 
  now also supported
* Included the ability to run LOHHLA on partial/sliced bam files, which only or at least 
  include the HLA region of chromosome 6, by allowing the user to input the number of 
  mapped reads for the tumor and normal `.bam` files. These two quantities are essential 
  to correctly compare the coverage between the two samples. I myself have used this 
  feature in combination with the TCGA bam slicing tool in order to estimate HLA LOH for 
  TCGA samples without having to download the entire `.bam` files for each patient.  For 
  this analysis, I estimated the total mapped reads using the file size of the `.bam` 
  file, which is directly accessible via TCGA .
* General usability, e.g. the ability to explicitly define the tumor and matched normal 
  `.bam` files.  The original currently expects as input a folder with exactly two `.bam` 
  files: the user defines the normal bam and the script infers that the other one must be 
  the tumor bam, I decided to to change this. The user can define both the normal and 
  tumor bams, allowing any kind of source file organization.
* More error checks and robustness -- although there's probably still room for improvement 
  here
* More informative tabular output. The added 'message' column will display why certain 
  alleles may have failed, so you don¿t necessarily have to go through the log files in 
  case a sample/allele fails (due to e.g. homozygosity or lack of coverage in the matched 
  normal bam)
* Solved some potential numerical problems in the code, e.g. I encountered a division 
  between two vectors of potentially unequal sizes, which lead to an error for some 
  samples
* Code readability and organization
* Some additional dependencies were included
* Plotting code has been minimally altered and remains untested in conjunction with the 
  rest of this codebase


A thank you goes out to Joris van der Haar for spotting some bugs I (Maarten Slagter) 
introduced.
Tested and developed in R 3.5 on Linux.

# README #

Immune evasion is a hallmark of cancer. Losing the ability to present productive tumor 
neoantigens could facilitate evasion from immune predation.  An integral part of 
neoantigen presentation is the HLA class I molecule, which presents epitopes to T-cells on 
the cell surface. Thus, loss of an HLA allele, resulting in HLA homozygosity, may be 
a mechanism of immune escape. However, the polymorphic nature of the HLA locus precludes 
accurate copy number calling using conventional copy number tools.  

Here, we present **LOHHLA**, **L**oss **O**f **H**eterozygosity in **H**uman **L**eukocyte 
**A**ntigen, a computational tool to evaluate HLA loss using next-generation sequencing 
data. 

### LICENCE ###

LOHHLA IS PROTECTED BY COPYRIGHT AND IS SUBJECT TO A PATENT APPLICATION.  THE TOOL IS 
PROVIDED “AS IS” FOR INTERNAL NON-COMMERCIAL ACADEMIC RESEARCH PURPOSES ONLY.   NO 
RESPONSIBILITY IS ACCEPTED FOR ANY LIABILITY ARISING FROM SUCH USE BY ANY THIRD PARTY.  
COMMERCIAL USE OF THIS TOOL FOR ANY PURPOSE IS NOT PERMITTED.  ALL COMMERCIAL USE OF THE 
TOOL INCLUDING TRANSFER TO A COMMERCIAL THIRD PARTY OR USE ON BEHALF OF A COMMERCIAL THIRD 
PARTY (INCLUDING BUT NOT LIMITED TO USE AS PART OF A SERVICE SUPPLIED TO ANY THIRD PARTY 
FOR FINANCIAL REWARD) REQUIRES A LICENSE.  FOR FURTHER INFORMATION PLEASE EMAIL Eileen 
Clark <eileen.clark@crick.ac.uk>.
 
 
### What do I need to install to run LOHHLA? ###

Please ensure a number of dependencies are first installed. These include:

* BEDTools (http://bedtools.readthedocs.io/en/latest/)
* SAMtools (http://samtools.sourceforge.net/)
* Novalign (http://www.novocraft.com/products/novoalign/)
* Picard (http://broadinstitute.github.io/picard/)
* R (https://www.r-project.org/about.html)

Within R, the following packages are required:

* seqinr (https://CRAN.R-project.org/package=seqinr)
* Biostrings (http://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* beeswarm (https://CRAN.R-project.org/package=beeswarm)
* zoo (https://cran.r-project.org/package=zoo)
* Rsamtools (http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
* dplyr (install.packages("dplyr"))
* naturalsort (install.packages("naturalsort"))

If not available locally, these packages will be attempted to be installed.

LOHHLA also requires an HLA fasta file. This can be obtained from Polysolver 
(http://archive.broadinstitute.org/cancer/cga/polysolver).

### How do I install LOHHLA? ###

To install LOHHLA, simply clone the repository:

git clone https://bitbucket.org/mcgranahanlab/lohhla.git

### How do I run LOHHLA? ###

LOHHLA is coded in R, and can be executed from the command line (Terminal, in 
Linux/UNIX/OSX, or Command Prompt in MS Windows) directly, or using a shell script (see 
example below).

USAGE: 

    Rscript /location/of/LOHHLA/script  [OPTIONS]

For a description of all the options, run: 

    Rscript /location/of/LOHHLA/script --help


### What is the output of LOHHLA? ###

LOHHLA produces multiple different files (see correct-example-out for an example). To 
determine HLA LOH in a given sample, the most relevant output is the file which ends 
'.HLAlossPrediction CI.xls'.  The most relavant columns are:

  HLA_A_type1  						   - the identity of allele 1
  HLA_A_type2  						   - the identity of allele 2
	Pval_unique  					     - this is a p-value relating to allelic imbalance 
	LossAllele      					 - this corresponds to the HLA allele that is subject to loss
	KeptAllele      					 - this corresponds to the HLA allele that is not subject to loss
	HLA_type1copyNum_withBAFBin          - the estimated raw copy number of HLA (allele 1)
	HLA_type2copyNum_withBAFBin          - the estimated raw copy number of HLA (allele 2)


For a full definition of the columns, see below, in each case whether the column should be used [use], or can be ignored [legacy]is indicated:

	region								 - the region or tumor sample [use]
	HLA_A_type1							 - the identity of allele 1 [use]
	HLA_A_type2							 - the identity of allele 2 [use]
	HLAtype1Log2MedianCoverage	         - the median LogR coverage across allele 1 [use] 
	HLAtype2Log2MedianCoverage	         - the median LogR coverage across allele 2 [use]
	HLAtype1Log2MedianCoverageAtSites	 - the median LogR coverage across allele 1, restricted to mismatch sites [use]
	HLAtype2Log2MedianCoverageAtSites	 - the median LogR coverage across allele 2, restricted to mismatch sites [use]
	HLA_type1copyNum_withoutBAF	         - estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withoutBAF_lower	 - lower 95% confidence interval of estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withoutBAF_upper	 - upper 95% confidence interval of estimated copy number of allele 1, without using BAF [legacy] 
	HLA_type1copyNum_withBAF	         - estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withBAF_lower	     - lower 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withBAF_upper	     - upper 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withoutBAF	         - estimated copy number of allele 2 without using BAF  [legacy] 
	HLA_type2copyNum_withoutBAF_lower	 - lower 95% confidence interval of estimated copy number of allele 2, without using BAF [legacy] 
	HLA_type2copyNum_withoutBAF_upper	 - upper 95% confidence interval of estimated copy number of allele 2, without using BAF [legacy] 
	HLA_type2copyNum_withBAF	         - estimated copy number of allele 2 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withBAF_lower	     - lower 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type2copyNum_withBAF_upper	     - upper 95% confidence interval of estimated copy number of allele 1 using BAF, without binning sites [legacy] 
	HLA_type1copyNum_withoutBAFBin	     - estimated copy number of allele 1 using binning, but without BAF [legacy]  
	HLA_type1copyNum_withoutBAFBin_lower - lower 95% confidence interval of estimated copy number of allele 1 using binning, but without BAF [legacy]  	
	HLA_type1copyNum_withoutBAFBin_upper - upper 95% confidence interval of estimated copy number of allele 1 using binning, but without BAF [legacy] 	
	HLA_type1copyNum_withBAFBin	         - estimated copy number of allele 1 using binning and BAF [use] 
	HLA_type1copyNum_withBAFBin_lower	 - lower 95% confidence interval of estimated copy number of allele 1 using binning and BAF [use] 
	HLA_type1copyNum_withBAFBin_upper	 - upper 95% confidence interval of estimated copy number of allele 1 using binning and BAF [use]  
	HLA_type2copyNum_withoutBAFBin	     - estimated copy number of allele 2 using binning, but without BAF [legacy]  
	HLA_type2copyNum_withoutBAFBin_lower - lower 95% confidence interval of estimated copy number of allele 2 using binning, but without BAF [legacy]	
	HLA_type2copyNum_withoutBAFBin_upper - upper 95% confidence interval of estimated copy number of allele 2 using BAF, without binning sites [legacy] 	
	HLA_type2copyNum_withBAFBin	         - estimated copy number of allele 2 using binning and BAF [use] 
	HLA_type2copyNum_withBAFBin_lower	 - lower 95% confidence interval of estimated copy number of allele 2 using binning and BAF [use] 
	HLA_type2copyNum_withBAFBin_upper	 - upper 95% confidence interval of estimated copy number of allele 2 using binning and BAF [use
	PVal                                 - p-value relating to difference in logR between allele 1 and allele 2 (paired t-test)[legacy]
	UnPairedPval	                     - p-value relating to difference in logR between allele 1 and allele 2 (unpaired t-test)[legacy]
	PVal_unique	                         - p-value relating to difference in logR between allele 1 and allele 2, ensuring each read only contributes once (paired t-test) [use]
	UnPairedPval_unique                  - p-value relating to difference in logR between allele 1 and allele 2, ensuring each read only contributes once (unpaired t-test) [use]
	LossAllele	                         - HLA allele that is present at lower frequency (potentially subject to loss) [use]
	KeptAllele                           - HLA allele that is present at higher frequency (potentially not subject to loss) [use]
	numMisMatchSitesCov                  - number of mismatch sites with sufficient coverage [use]
	propSupportiveSites                  - proportion of missmatch sites that are consistent with loss or allelic imbalance [use]


### How can I test if LOHHLA is working? ###

Example data is included in the LOHHLA repository. To run LOHHLA on the example dataset, 
alter the "example.sh" script to match your local file structure and ensure the requisite 
dependencies are available / loaded.  The `--HLAfastaLoc`, `--gatkDir`, and `--novoDir` 
file paths should also be updated to the corresponding locations.  File paths must be full 
paths. Run "example.sh" and the output should match that found in the 
`correct-example-out` directory provided.  All BAM files (normal and tumour) should be 
found in or linked to the same directory.

### Who do I talk to? ###

If you have any issues with LOHHLA, please send an email to lohhla@gmail.com

### How do I cite LOHHLA ? ###

If you use LOHHLA in your research, please cite the following paper:

McGranahan et al., Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell (2017), https://doi.org/10.1016/j.cell.2017.10.001

