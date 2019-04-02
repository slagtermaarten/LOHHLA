library(optparse)

option_list = list(
  make_option(c("-id", "--patientId"), type="character", default=NULL,
              help="patient ID", metavar="character"),
  make_option(c("-o", "--outputDir"), type="character", default=NULL,
              help="location of output directory", metavar="character"),
  make_option(c("-tBAM", "--tumorBAMfile"), type="character", default=NULL,
              help="tumor BAM file", metavar="character"),
  make_option(c("-nBAM", "--normalBAMfile"), type="character", default=NULL,
              help="normal BAM file\n\t\tcan be FALSE to run without normal sample", metavar="character"),
  make_option(c("-BAM", "--BAMDir"), type="character", default=NULL,
              help="location of all BAMs to test", metavar="character"),
  make_option(c("-hla", "--hlaPath"), type="character", default=NULL,
              help="location to patient HLA calls", metavar="character"),
  make_option(c("-hlaLoc", "--HLAfastaLoc"), type="character", 
                default='~/lohhla/data/hla_all.fasta',
                help="location of HLA FASTA [default= %default]", 
                metavar="character"),
  make_option(c("-LOHHLA_loc"), type="character", 
                default='~/libs/LOHHLA',
                help="location of LOHHLA R library [default= %default]", 
                metavar="character"),
  make_option(c("-cn", "--CopyNumLoc"), type="character", default="FALSE",
              help="location to patient purity and ploidy output\n\t\tcan be FALSE to only estimate allelic imbalance", metavar="character"),
  make_option(c("-ov", "--overrideDir"), type="character", default='FALSE',
              help="location of flagstat information if already run [default= %default]", metavar="character"),
  make_option(c("-mc", "--minCoverageFilter"), type="numeric", default=30,
              help="minimum coverage at mismatch site [default= %default]", metavar="character"),
  make_option(c("-kmer", "--kmerSize"), type="numeric", default=50,
              help="size of kmers to fish with [default= %default]", metavar="character"),
  make_option(c("-mm", "--numMisMatch"), type="numeric", default=1,
              help="number of mismatches allowed in read to map to HLA allele [default= %default]", metavar="character"),
  make_option(c("-ms", "--mappingStep"), type="logical", default=TRUE,
              help="does mapping to HLA alleles need to be done [default= %default]", metavar="character"),
  make_option(c("-fs", "--fishingStep"), type="logical", default=TRUE,
              help="if mapping is performed, also look for fished reads matching kmers of size kmerSize [default= %default]", metavar="character"),
  make_option(c("-ps", "--plottingStep"), type="logical", default=TRUE,
              help="are plots made [default= %default]", metavar="character"),
  make_option(c("-cs", "--coverageStep"), type="logical", default=TRUE,
              help="are coverage differences analyzed [default= %default]", metavar="character"),
  make_option(c("-cu", "--cleanUp"), type="logical", default=TRUE,
              help="remove temporary files [default= %default]", metavar="character"),
  make_option(c("-no", "--novoDir"), type="character", default='',
              help="path to novoalign executable [default= %default]", metavar="character"),
  make_option(c("-na", "--novoThreads"), type="integer", default=min(as.integer(system('nproc', intern = T)), 8),
              help="amount of threads to be used by novoalign [default= %default]", metavar="character"),
  make_option(c("-ga", "--gatkDir"), type="character", default='',
              help="path to GATK executable [default= %default]", metavar="character"),
  make_option(c("-ex", "--HLAexonLoc"), type="character", default='~/lohhla/data/hla.dat',
              help="HLA exon boundaries for plotting [default= %default]", metavar="character"),
  make_option(c("-w", "--ignoreWarnings"), type="logical", default=TRUE,
              help="continue running with warnings [default= %default]", metavar="character"),
  make_option(c("-as", "--genomeAssembly"), type="character", default='grch38',
              help="specify genome assembly (hg19 or grch38) [default= %default]", metavar="character"),
  make_option(c("-jl", "--jellyFish"), type="character", default='jellyfish',
              help="specify location of jellyfish binary", metavar="character"),
  make_option(c("-bt", "--bedtools"), type="character", default='bedtools',
              help="specify location of bedtools binary", metavar="character"),
  make_option(c("-st", "--samtools"), type="character", default='samtools',
              help="specify location of samtools binary", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## TODO infer this automatically or just configure this to a package
if (!dir.exists(opt$LOHHLA_loc)) {
  stop('LOHHLA_loc command line argument seems to be incorrect')
}
source(file.path(opt$LOHHLA_loc, 'LOHHLA_funcs.R'))

# strip trailing / from all parameters in opt
for (i in 1:length(opt)) {
  if (substr(opt[i], nchar(opt[i]), nchar(opt[i])) == '/') {
    opt[i] <- substr(opt[i], 1, nchar(opt[i]) - 1)
  }
}
print(opt)

run_LOHHLA(opt)
