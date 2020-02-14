if (!require('yaml')) {
  install.packages('yaml')
}
library(yaml)
## Read paths from yaml so they won't have to be included in git repo
config <- yaml.load_file(input = '~/libs/LOHHLA/config.yaml')
source(file.path(config$LOHHLA_dir, 'functions.R'))
example_root <- file.path(config$LOHHLA_dir, 'example-file')
stored_res_dir <- file.path(example_root, 'correct-example-out')
new_res_dir <- file.path(example_root, 'out')
use_original_mapping <- F

if (use_original_mapping) {
  ## Don't redo mapping steps in order to investigate whether downstream
  ## analysis steps could be causing a discrepancy with original output
  new_res_om <- file.path(example_root, 'out_old_mapping')
  # dir.create(new_res_om, showWarnings = F)
  # list.files(new_res_om)
  # list.files(stored_res_dir)
  system(sprintf('rsync -av %s/ %s', stored_res_dir, new_res_om))
  system('rm %s', 
         file.path(new_res_om, 'example.10.DNA.HLAlossPrediction_CI.xls'))
  new_wd <- new_res_om
} else {
  new_wd <- stored_res_dir
}

opt <- list(
  patientId = 'example',
  performAlignment = F,
  outputDir = new_wd,
  normalBAMfile = file.path(example_root, 'bam', 'example_BS_GL_sorted.bam'),
  tumorBAMfile = file.path(example_root, 'bam', 'example_tumor_sorted.bam'),
  BAMDir = file.path(example_root, 'bam'),
  hlaPath = file.path(example_root, 'hlas'),
  HLAfastaLoc = file.path(config$LOHHLA_dir, 'data', 'example.patient.hlaFasta.fa'),
  HLAexonLoc = file.path(config$LOHHLA_dir, 'data', 'hla.dat'),
  CopyNumLoc = file.path(config$LOHHLA_dir, 'example-test', 'solutions.txt'),
  genomeAssembly = 'hg19',
  mappingStep = !use_original_mapping %||% T,
  coverageStep = !use_original_mapping %||% T,
  plottingStep = FALSE,
  minCoverageFilter = 10,
  debug = T,
  forceRedo = !use_original_mapping %||% T,
  requirePairedReads = F,
  fishingStep = TRUE,
  cleanUp = FALSE,
  jellyFish = config$jellyFish,
  gatkDir = config$gatkDir,
  novoDir = config$novoDir)

if (F) {
  ## REDUNDANT
  ## Copy over old results to see whether their difference causes final results
  ## to be different
  file_id <-
    c('example_tumor_sorted/example_tumor_sorted.chr6region.patient.reference.hlas.sam',
      'example_BS_GL_sorted/example_BS_GL_sorted.chr6region.patient.reference.hlas.sam')
  for (f in file_id) {
    file.copy(file.path(stored_res_dir, f), file.path(new_res_dir, f))
  }
}

run_LOHHLA(opt)

if (T) {
  latest_log <- list.files(new_wd, pattern = 'running', 
                           full.names = T) %>%
    .[length(.)]
  # maartenutils::less(latest_log)
  ## Check results against expected results
  new_fn <- list.files(new_wd, pattern = 'HLAlossPrediction_CI', 
                       full.names = T) %>%
    .[length(.)]
  print(file.mtime(new_fn))
  new_res <- new_fn %>%
    read.table(header = T) %>%
    .[!is.na(.$HLA_A_type1), ] %>%
    { tryCatch(., error = function(e) { print(e) }) }

  original_res <-
    read.table(file.path(stored_res_dir, 
                         'example.10.DNA.HLAlossPrediction_CI.xls'),
               header = T)

  shared_cols <- intersect(colnames(original_res), colnames(new_res))
  original_res <- original_res[, shared_cols, drop = F]
  new_res <- new_res[, shared_cols, drop = F]
  combined <- rbind(original_res, new_res, fill = F)
  t_combined <- combined[!is.na(combined$HLA_A_type1), ] %>% t
  t_combined <- cbind(as.data.frame(t_combined), 
                      'congruence' = apply(t_combined, 1, 
                                           function(x) length(unique(x)) == 1))
  print(t_combined)
}
