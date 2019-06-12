## {{{ Install required libraries if required
for (p in c('seqinr', 'Biostrings', 'beeswarm', 'zoo', 'Rsamtools', 'dplyr',
    'naturalsort', 'glue', 'magrittr', 'tools', 'purrr')) {
  if (!require(p, character.only = T)) {
    install.packages(p, character = T)
  }
  require(p, quietly = T, character.only = T)
}
## }}}

## {{{ Helper functions
#' Write a tsv to disk
#'
#'
write_tsv <- function(object, filename, append = F,
  col.names = T, row.names = F) {
  if (append == T && !file.exists(filename)) {
    dir.create(dirname(filename), recursive = T, showWarnings = F)
    system(glue('touch {filename}'), intern = T)
  }
  # less(filename)
  tryCatch(
    write.table(object, file = filename, sep = '\t',
      quote = FALSE, row.names = row.names,
      append = append, col.names = col.names)
    , error = function(e) { print(e); browser() },
    warning = function(e) { })
}


#' Replace NA elements of a vector with a replacement
#'
#'
repl_NA <- function(vec, replacement = NA) {
  sapply(vec, function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x))
      return(replacement)
    else
      return(x)
  })
}


#' Remove NA elements from a vector
#'
#'
remove_NA <- function(vec) {
  vec[!is.na(vec)]
}


#' Column bind a set of potentially unequal length vectors
#'
#'
cbind_uneven <- function(...) {
  max_length <- max(sapply(list(...), length), na.rm = T)
  extended_items <- lapply(list(...), function(x) {
    N_missing_elems <- max_length - length(x)
    c(x, rep(NA, N_missing_elems))
  })
  do.call(cbind, extended_items)
}


#' Divide two vectors and log2 transform the result, only using shared elements
#'
#'
element_divide_vector <- function(
  a, b, sites = NULL, mult_factor = 1, extra = 0.0001) {
  stopifnot(is.vector(a))
  stopifnot(is.vector(b))
  stopifnot(names(a) != NULL)
  stopifnot(names(b) != NULL)

  if (length(a) == 0 || length(b) == 0) {
    return(NA)
  }

  if (is.null(sites)) {
    sites <- intersect(names(a), names(b))
  }
  shared_sites <- sites %>%
    as.character %>%
    intersect(names(a)) %>%
    intersect(names(b))

  res <- tryCatch(data.frame(
    'pos' = as.numeric(shared_sites),
    'logR' = log2(a[shared_sites] / b[shared_sites]) + extra * mult_factor,
    row.names = NULL
  ), error = function(e) { browser() }, warning = function(e) { browser() })

  return(res)
}


#' Harmonize the colnames of two matrix objects such that they can be
#' row bound
#'
#'
harmonize_colnames <- function(mat1, mat2) {
  if (!is.null(mat1)) {
    missing_cols <- setdiff(colnames(mat1), colnames(mat2))
    if (length(missing_cols) > 0) {
      for (coln in missing_cols) {
        mat2 <- cbind(mat2, coln = NA)
        colnames(mat2)[length(colnames(mat2))] <- coln
      }
    }

    new_cols <- setdiff(colnames(mat2), colnames(mat1))
    if (length(new_cols) > 0) {
      for (coln in new_cols) {
        mat1 <- cbind(mat1, coln = NA)
        colnames(mat1)[length(colnames(mat1))] <- coln
      }
    }

    if (!all(colnames(mat2) == colnames(mat1))) {
      ## Order identically
      idx <- match(colnames(mat2), colnames(mat1))
      mat2 <- mat2[, idx]
    }
  }
  return(list(mat1, mat2))
}


document_params <- function(params, log_fn) {
  msg <- "\n#######################\ \n####### Inputs ######## \n#######################\n"

  for (i in seq_along(params)) {
    tmp <- c()
    for (j in 1:length(params[[i]])) {
      tmp <- paste(tmp, params[[i]][j], "\t", sep = "")
    }
    msg <- paste(msg, names(params)[i], "\t", tmp, "\n", sep = "")
  }

  write_tsv(msg, log_fn, append = TRUE)
}


PasteVector <- function(v, sep = "") {
  vt <- v[1]
  if (length(v) > 1) {
    for (g in 2:length(v)) {
      vt <- paste(vt, v[g], sep = sep)
    }
  }
  vt <- paste(vt, ' EnD', sep = '')
  out.v <- sub(' EnD', '', vt)
  out.v <- sub('NA , ', '', out.v)
  out.v <- sub(' , NA', '', out.v)
  out.v <- sub(' , NA , ', ' , ', out.v)
  return(out.v)
}


create_logger <- function(log_fn) {
  ## See Hadley Wickham's Advanced-R on functionals to see how this
  ## works if unclear
  force(log_fn)
  if (!file.exists(log_fn)) {
    dir.create(dirname(log_fn), recursive = T, showWarnings = F)
    system(sprintf('touch %s', log_fn), intern = T)
  }
  logger <- function(...) {
    msg <- paste(..., collapse = ' ')
    f_msg <- sprintf('%s - %s\n', date(), msg)
    message(f_msg)
    write.table(f_msg, file = log_fn,
      quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  return(logger)
}
# logger <- create_logger('test.log')
# logger('aa', 'bb')


identify_kmers <- function(workDir, kmerSize, HLAfastaLoc, jellyFish) {
  setwd(workDir)
  cmd <- paste(jellyFish, 'count -m', kmerSize, '-s 100M -t 5', HLAfastaLoc)
  system(cmd)
  ## mer_counts_dumps.fa is a fasta formatted file like the following:
  ## >1
  ## ACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACGCCGAGCTTGTGGA
  ## >1372
  ## CCTCCTTTCCCAGAGCCATCTTCCCAGTCCACCATCCCCATCGTGGGCAT
  ## >1
  ## GGGCGCCGTGGATGGAGCAGGAGGGGCCGGAGTATTGGGACCGGAACACA
  ## In which the >\\d indicate k-mer counts and the DNA sequences are the
  ## 50-mers as identified using jellyfish
  cmd <- paste(jellyFish, 'dump mer_counts.jf > mer_counts_dumps.fa')
  system(cmd)

  ## Strip away the k-mer counts to only retain the unique k-mers
  cmd <- paste('grep --invert-match \\> mer_counts_dumps.fa > ',
    kmerSize, 'mer_uniq', sep = '')
  system(cmd)
}


get.partially.matching.reads <- function(workDir, sample_dir, BAMDir, BAMfile,
  kmerSize, gatkDir) {
  kmerFile <- paste(workDir, '/', kmerSize, 'mer_uniq', sep = '')

  # add header
  cmd <- paste(samtools, ' view -H ', BAMDir, '/', BAMfile, ' > ',
               sample_dir, '/fished.sam', sep = '')
  system(cmd)

  # fish partially matching reads
  cmd <- paste(samtools, ' view ', BAMDir, '/', BAMfile, ' | grep -F -f ',
               kmerFile, ' >> ', sample_dir, '/fished.sam', sep = '')
  system(cmd)

  ## convert to fastq
  cmd <- paste('java -jar ', gatkDir, '/SamToFastq.jar',
    ' I=', sample_dir, '/fished.sam',
    ' F=', sample_dir, '/fished.1.fastq',
    ' F2=', sample_dir, '/fished.2.fastq VALIDATION_STRINGENCY=SILENT', sep = '')
  system(cmd)
}


combine.fastqs <- function(chr6.f1, chr6.f2, fished.f1, fished.f2) {
  chr6.f1.seq <- read.table(chr6.f1, sep = '\t', stringsAsFactors = FALSE, comment.char = '', quote = '')
  chr6.f2.seq <- read.table(chr6.f2, sep = '\t', stringsAsFactors = FALSE, comment.char = '', quote = '')
  fished.f1.seq <- tryCatch(read.table(fished.f1, sep = '\t',
    stringsAsFactors = FALSE, comment.char = '', quote = ''),
    error = function(e) { print(e); NULL })
  fished.f2.seq <- tryCatch(read.table(fished.f2, sep = '\t',
      stringsAsFactors = FALSE, comment.char = '', quote = ''),
    error = function(e) { print(e); NULL })
  fished.f1.names <- fished.f1.seq[seq(1, nrow(fished.f1.seq), by = 4), ]
  fished.f2.names <- fished.f2.seq[seq(1, nrow(fished.f2.seq), by = 4), ]
  fished.f1.names.nodup <- fished.f1.names[-which(fished.f1.names %in% chr6.f1.seq$V1)]
  fished.f2.names.nodup <- fished.f2.names[-which(fished.f2.names %in% chr6.f2.seq$V1)]
  fished.f1.seq.toadd <- lapply(which(fished.f1.seq$V1 %in% fished.f1.names.nodup), FUN = function(x) {
    return(fished.f1.seq[x:(x+3), ])
    })
  fished.f1.seq.toadd <- unlist(fished.f1.seq.toadd)
  fished.f2.seq.toadd <- lapply(which(fished.f2.seq$V1 %in% fished.f2.names.nodup), FUN = function(x) {
    return(fished.f2.seq[x:(x+3), ])
    })
  fished.f2.seq.toadd <- unlist(fished.f2.seq.toadd)

  new.chr6.f1.seq <- c(chr6.f1.seq$V1, fished.f1.seq.toadd)
  new.chr6.f2.seq <- c(chr6.f2.seq$V1, fished.f2.seq.toadd)

  write_tsv(new.chr6.f1.seq, chr6.f1, row.names = FALSE, col.names = FALSE)
  write_tsv(new.chr6.f2.seq, chr6.f2, row.names = FALSE, col.names = FALSE)
}


count_events <- function(BAMfile, n) {
  x <- scanBam(BAMfile, index = BAMfile,
    param = ScanBamParam(what = scanBamWhat(), tag = 'NM'))
  readIDs <- x[[1]][['qname']]
  if (length(readIDs) == 0) return(NULL)
  cigar <- x[[1]][['cigar']]
  editDistance <- unlist(x[[1]][['tag']])
  insertionCount <- sapply(cigar, FUN = function(boop) {
    length(grep(pattern = 'I', x = unlist(strsplit(boop, split = ''))))
  })
  deletionCount <- sapply(cigar, FUN = function(boop) {
    length(grep(pattern = 'D', x = unlist(strsplit(boop, split = ''))))
  })
  indelTotals <- sapply(cigar, FUN = function(boop) {
    tmp <- unlist(strsplit( gsub('([0-9]+)', '~\\1~', boop), '~' ))
    Is <- grep(pattern = 'I', x = tmp)
    Ds <- grep(pattern = 'D', x = tmp)
    total <- sum(as.numeric(tmp[(Is-1)])) + sum(as.numeric(tmp[Ds-1]))
    return(total)
  })
  misMatchCount <- tryCatch(editDistance - indelTotals,
    error = function(e) { print(e); browser() })
  eventCount <- misMatchCount + insertionCount + deletionCount
  names(eventCount) <- 1:length(eventCount)
  passed <- eventCount[which(eventCount <= n)]
  y <- readIDs[as.numeric(names(passed))]
  y <- names(table(y)[which(table(y) == 2)])
  return(y)
}


dont.count.twice <- function(BAMfile1, BAMfile2,
  normalBAMfile1, normalBAMfile2) {

  x <- scanBam(BAMfile1)
  readIDs_x <- x[[1]][['qname']]
  y <- scanBam(BAMfile2)
  readIDs_y <- y[[1]][['qname']]

  table(readIDs_x %in% readIDs_y)

  xn <- scanBam(normalBAMfile1)
  readIDs_xn <- xn[[1]][['qname']]
  yn <- scanBam(normalBAMfile2)
  readIDs_yn <- yn[[1]][['qname']]

  remove.from.tumor <- readIDs_x[which(readIDs_x %in% readIDs_y)]
  remove.from.normal <- readIDs_xn[which(readIDs_xn %in% readIDs_yn)]

  out <- list(remove.from.tumor, remove.from.normal)
  names(out) <- c('remove.from.tumor', 'remove.from.normal')
  return(out)
}


getMisMatchPositionsPairwiseAlignment <- function(alignment,
  chunksize = 60, returnlist = FALSE) {

  ## Aligned residues
  seq1_ar <- unlist(strsplit(as.character(pattern(alignment)), split = ''))
  seq2_ar <- unlist(strsplit(as.character(subject(alignment)), split = ''))

  k <- 1
  seq1Positions <- c()
  for (char in seq1_ar) {
    if (char %in% c('C', 'G', 'A', 'T')) {
      seq1Positions <- c(seq1Positions, k)
      k <- k + 1
      next;
    }
    if (char %in% c('-')) {
      seq1Positions <- c(seq1Positions, k)
      next;
    }
  }

  k <- 1
  seq2Positions <- c()
  for (char in seq2_ar) {
    if (char %in% c('C', 'G', 'A', 'T')) {
      seq2Positions <- c(seq2Positions, k)
      k <- k + 1
      next;
    }
    if (char %in% c('-')) {
      seq2Positions <- c(seq2Positions, k)
      next;
    }
  }

  mm_bools <- seq1_ar != seq2_ar
  diffSeq1 <- seq1Positions[mm_bools]
  diffSeq2 <- seq2Positions[mm_bools]

  diffType1 <- rep(1, length(diffSeq1))
  diffType1[which(seq1_ar[mm_bools] %in% '-')] <- 2

  diffType2 <- rep(1, length(diffSeq2))
  diffType2[which(seq2_ar[mm_bools] %in% '-')] <- 2

  out <- list(
    diffSeq1 = diffSeq1,
    diffSeq2 = diffSeq2,
    diffType1 = diffType1,
    diffType2 = diffType2
  )
  return(out)
}


getUniqMapReads <- function(workDir, BAMDir, override = FALSE,
  tumorBAMfile, normalBAMfile, overrideDir = NULL) {
  if (!override) {
    outDir <- paste(workDir, '/flagstat/', sep = '')
    if (!file.exists(outDir)) {
      if (!dir.create(outDir, recursive = TRUE)) {
        msg <- paste("Unable to create directory: ", outDir, "!\n", sep = '')
        logger(msg)
        stop(msg)
      }
    }

    BAMs <- list.files(BAMDir, pattern = 'bam$', full.names = TRUE)

    for (BAM in c(tumorBAMfile, normalBAMfile)) {
      sample <- unlist(strsplit(BAM, split = '/')) %>% { .[length(.)] }
      flagstat <- paste0(outDir, sample, '.proc.flagstat')
      cmd <- paste0(samtools, ' flagstat ', BAM, ' > ', flagstat)
      if (file.exists(flagstat)) {
        logger(flagstat, 'already exists, skipping')
      } else {
        system(cmd)
      }
    }
  } else {
    outDir <- overrideDir
  }

  flagStatRegions <- list.files(outDir, pattern = ".proc.flagstat$") %>%
    { grep('filtered', ., invert = T, value = T) }

  if (length(flagStatRegions) == 0) {
    stop('Either run flagstat or do not override.')
  }

  UniqMapReads <- lapply(inverseNames(baseName(flagStatRegions)),
    function(x) {
      fn <- paste0(outDir, x)
      as.numeric(read.table(fn, stringsAsFactors = FALSE, header = FALSE,
                            nrows = 1)[, 1])
    })
  return(UniqMapReads)
}


funCalcN_withBAF <- function(logR, BAF, psi, rho, gamma = 1) {
  avg_ploidy <- (1 - rho) * 2 + rho * psi
  nA <- (rho - 1 + BAF * 2^(logR/gamma) * avg_ploidy) / rho
  nB <- (rho - 1 - (BAF - 1) * 2^(logR/gamma) * avg_ploidy) / rho
  return(cbind(nA, nB))
}


funCalcN_withoutBAF <- function(logR, psi, rho, gamma) {
  ((((1 - rho) + rho * psi / 2)) * 2^(logR/gamma) - (1 - rho)) / rho
}


t.test.NA <- function(...) {
  tryCatch(t.test(...), error = function(e)
    list('stat' = NA, 'conf.int' = c(NA, NA)))
}


#' Get exon positions for HLA
#'
#'
getHlaExons <- function(allele, pickAllele = 'all') {
  alleleExonFormatted <- strsplit(allele, split = "_") %>%
    unlist %>%
    { paste('HLA-', .[2], '\\*',
            paste(.[3:length(.)], collapse = ':'), sep = '') } %>%
    toupper()
  ## Retrieve first matching allele
  cmd <- paste("awk \'/^DE   ", alleleExonFormatted,
               "/ {p = 1}; p; /^SQ   Sequence/ {p = 0}\' ", HLAexonLoc, sep  = '')
  awkResult <- system(cmd, intern = TRUE)
  if (length(awkResult) == 0) return(NULL)

  if (pickAllele == 'first') {
    awkResult %<>% { .[1:(which(grepl('^SQ', .))[1])] }
  } else if (pickAllele == 'largest') {
    ## TODO implement me
    # awkResult %>% { which.max(diff(which(grepl('^SQ', .)))) + 1 }
  } else if (pickAllele == 'all') {
    ## No selection
  }

  hlaExons <- grep('^FT   exon', awkResult, value = TRUE) %>%
    { strsplit(., split = " ") } %>%
    { do.call(rbind, .) } %>%
    { .[, ncol(.)] } %>%
    unique %>%
    naturalsort

  if (length(hlaExons) != 0) {
    exonTable <- c()
    for (i in 1:length(hlaExons)) {
      exonTable <- rbind(exonTable,
        unlist(strsplit(hlaExons[i], split = '\\.\\.')))
    }
  }
  return(list(hlaExons, exonTable))
}


#' Get base name of file and reduce it to everything leading up to first dot
#'
#'
baseName <- function(fns) {
  sapply(fns, function(fn) unlist(strsplit(basename(fn), split = "\\."))[1])
}


#' Exchange a character vector's names and values
#'
#'
inverseNames <- function(vec) {
  stopifnot(is.character(vec))
  stopifnot(is.character(names(vec)))
  setNames(names(vec), vec)
}
## }}}

run_LOHHLA <- function(opt) {
  ## {{{ Preamble
  ## Strip trailing slashes away from directories
  opt <- lapply(opt, function(x) if (is.character(x)) gsub('/$', '', x) else x)

  # Check if variables are already defined in current environment.
  # Replace undefined ones with command line arguments or defaults
  patientId <- opt$patientId
  workDir <- outputDir <- opt$outputDir
  normalBAMfile <- opt$normalBAMfile
  tumorBAMfile <- opt$tumorBAMfile
  error_msg <- ''
  # example of invalid file name (as in LOHHLA code cannot handle this currently)
  # grepl('^.+\\..+\\.bam$', "normal_hla_loci.sorted.bam")
  # test for valid file name
  ## Some input checking
  if (grepl('^.+\\..+\\.bam$', basename(normalBAMfile)) ||
      grepl('^.+\\..+\\.bam$', basename(tumorBAMfile))) {
    stop('LOHHLA assumes absence of dots in bam filenames, please change names')
  }
  # if (grepl('hla|HLA', basename(normalBAMfile))) {
  #   warning('LOHHLA cannot handle the word "hla" filenames, proceed at own risk')
  # }

  if (!is.null(opt$HLAfastaLoc) && !file.exists(opt$HLAfastaLoc)) {
    opt$HLAfastaLoc <- NULL
  }

  ## Unpack command line arguments and check the paths
  BAMDir <- opt$BAMDir
  stopifnot(dir.exists(BAMDir))
  hlaPath <- opt$hlaPath
  stopifnot(file.exists(hlaPath))
  HLAfastaLoc <- opt$HLAfastaLoc
  stopifnot(file.exists(HLAfastaLoc))
  CopyNumLoc <- opt$CopyNumLoc
  stopifnot(file.exists(CopyNumLoc))
  overrideDir <- opt$overrideDir
  minCoverageFilter <- opt$minCoverageFilter
  numMisMatch <- opt$numMisMatch
  mappingStep <- opt$mappingStep
  cleanUp <- opt$cleanUp
  gatkDir <- opt$gatkDir
  stopifnot(dir.exists(gatkDir))
  kmerSize <- opt$kmerSize
  fishingStep <- opt$fishingStep
  coverageStep <- opt$coverageStep
  plottingStep <- opt$plottingStep
  ignoreWarnings <- opt$ignoreWarnings
  HLAexonLoc <- opt$HLAexonLoc
  stopifnot(file.exists(HLAexonLoc))
  novoDir <- opt$novoDir
  stopifnot(dir.exists(novoDir))
  stopifnot('novoalign' %in% list.files(novoDir, all.files = T))
  stopifnot('novoindex' %in% list.files(novoDir, all.files = T))
  novoThreads <- opt$novoThreads
  genomeAssembly <- tolower(opt$genomeAssembly)
  jellyFish <- opt$jellyFish
  stopifnot(file.exists(jellyFish))
  bedtools <- opt$bedtools
  stopifnot(file.exists(bedtools))
  samtools <- opt$samtools
  stopifnot(file.exists(samtools))
  fnExt <- opt$fnExt

  if (is.null(BAMDir) || is.null(outputDir) ||
      is.null(normalBAMfile) ||
      is.null(tumorBAMfile) ||
      is.null(hlaPath) || is.null(HLAfastaLoc)) {
    if (interactive()) print_help(opt_parser)
    stop('Missing arguments.\n', call. = FALSE)
  }

  figureDir <- file.path(outputDir, 'Figures')
  log_fn <- paste(outputDir, '/running.hla.loh.exome@',
    gsub(':', '-', gsub(' +', '_', date())), '.log', sep = '')
  logger <<- create_logger(log_fn)

  runWithNormal <- any(normalBAMfile != 'FALSE')
  extractNONmismatchReads <- TRUE
  extractUniqueReads <- TRUE
  performIntegerCopyNum <- TRUE
  useLogRbin <- TRUE
  if (CopyNumLoc == 'FALSE') {
    performIntegerCopyNum <- FALSE
    useLogRbin <- FALSE
  }
  override <- ifelse(overrideDir == 'FALSE', yes = FALSE, no = TRUE)
  gamma <- 1
  binSize <- 150

  if (ignoreWarnings) {
    howToWarn <- function(x) warning(x, call. = F, immediate. = T)
  } else {
    howToWarn <- function(x) stop(x)
  }
  howToWarn <- function(x) { }

  if (!dir.exists(workDir)) {
    dir.create(workDir, recursive = TRUE)
  }

  if (!dir.exists(figureDir)) {
    dir.create(figureDir, recursive = TRUE)
  }

  params <- list(patientId, workDir, hlaPath, normalBAMfile, tumorBAMfile,
    BAMDir, HLAfastaLoc, CopyNumLoc, gatkDir, novoDir)
  names(params) <- c('patientId', 'workDir', 'hlaPath', 'normalBAMfile',
    'tumorBAMfile', 'BAMDir', 'HLAfastaLoc', 'CopyNumLoc', 'gatkDir', 'novoDir')
  document_params(params, log_fn)

  BAMfiles <- c(basename(tumorBAMfile), basename(normalBAMfile))
  samples <- sapply(BAMfiles,
    function(x) unlist(strsplit(x, split = '.bam'))[1])
  normalName <- samples[basename(normalBAMfile)]
  tumorName <- samples[basename(tumorBAMfile)]
  ## Preamble }}}

  ## {{{ Read HLA alleles
  hlaAlleles <- read.table(hlaPath, sep = '\t', header = FALSE, as.is = TRUE)
  if (ncol(hlaAlleles) == 3) {
    hlaAlleles <- unique(sort(c(hlaAlleles$V2, hlaAlleles$V3)))
  } else if (ncol(hlaAlleles) == 1) {
    hlaAlleles <- unique(sort(hlaAlleles$V1))
  } else {
    stop('Unexpected parsing of hlaPath')
  }

  hlaFasta <- read.fasta(HLAfastaLoc)
  if (!all(hlaAlleles %in% names(hlaFasta))) {
    missing_hlas <- setdiff(hlaAlleles, names(hlaFasta))
    msg <- paste('Missing HLAs from FASTA: ',
      paste(missing_hlas, collapse = ', '),
      '. Trying to find alternative(s).', sep = '')
    howToWarn(msg)
    logger(msg)

    ## Attempt to subsitute HLAs not found in ref with more specific HLAs
    for (hla in missing_hlas) {
      alt <- grep(pattern = hla, x = names(hlaFasta), value = TRUE)[1]
      if (!is.na(alt)) {
        msg <- paste('Replacing: ', hla, ' with ', alt,  '.', sep = '')
        howToWarn(msg)
        logger(msg)
        hlaAlleles[which(hlaAlleles == hla)] <- alt
      } else {
        msg <- glue::glue('Could not find replacement for {hla}, \\
          removing it from HLA repertoire')
        howToWarn(msg)
        logger(msg)
        hlaAlleles <- setdiff(hlaAlleles, hla)
      }
    }
  }

  ## check for homozygous alleles here to save time on mapping step.
  ## also figure out if hla names will be uniformly 'hla_x'
  if (length(grep('hla_a', x = hlaAlleles)) == 1) {
    msg <- 'Homozygous for HLA-A -- not going to see any LOH here'
    logger(msg)
    howToWarn(msg)
    hlaAlleles <- hlaAlleles[-grep('hla_a', x = hlaAlleles)]
  }

  if (length(grep('hla_b', x = hlaAlleles)) == 1) {
    msg <- 'Homozygous for HLA-B -- not going to see any LOH here'
    logger(msg)
    howToWarn(msg)
    hlaAlleles <- hlaAlleles[-grep('hla_b', x = hlaAlleles)]
  }

  if (length(grep('hla_c', x = hlaAlleles)) == 1) {
    msg <- 'Homozygous for HLA-C -- not going to see any LOH here'
    logger(msg)
    howToWarn(msg)
    hlaAlleles <- hlaAlleles[-grep('hla_c', x = hlaAlleles)]
  }

  if (length(hlaAlleles) == 0) {
    ## Turn mapping of because there's nothing to map anyway
    mappingStep <- F
    coverageStep <- F
  } 
  ### }}}

  ## {{{ Pre-mapping step: double check whether we actually want to do the
  ## mappingStep or that we can skip it.
  if (mappingStep) {
    all_exist <- TRUE
    for (BAMfile in BAMfiles) {
      BAMid <- unlist(strsplit(BAMfile, split = '.bam'))[1]
      sample_dir <- paste(workDir, '/', BAMid, sep = '')
      for (allele in hlaAlleles) {
        bai <- paste0(sample_dir, '/', BAMid, '.type.', allele,
          '.filtered.bam.bai')
        all_exist <- all_exist && file.exists(bai)
      }
      ## We don't need to keep on checking if it's already FALSE
      if (!all_exist) break
    }
    mappingStep <- !all_exist
    if (!mappingStep) {
      msg <- glue::glue('Observed all allele .filtered.bam.bai files, \\
        so I assume the mapping step is already done! SKIPPING!')
      logger(msg)
    }
  }
  ## }}}

  ## {{{ Mapping step
  if (mappingStep) {
    ## generate patient reference fasta
    logger('generate patient reference fasta')
    patient.hlaFasta <- hlaFasta[hlaAlleles]
    patient_hla_fasta_fn <-
      paste(workDir, '/', patientId, '.patient.hlaFasta.fa', sep = '')
    patient_hla_fasta_index_fn <-
      paste(workDir, '/', patientId, '.patient.hlaFasta.nix', sep = '')
    write.fasta(patient.hlaFasta, file = patient_hla_fasta_fn,
      names = names(patient.hlaFasta))

    ## nix file for patient reference fasta -- novoalign
    logger('generating nix file for patient reference fasta')
    novoindexCMD <- paste0(novoDir, '/novoindex -t ', novoThreads, ' ',
      patient_hla_fasta_index_fn, ' ', patient_hla_fasta_fn)
    logger(novoindexCMD)
    system(novoindexCMD)

    if (fishingStep) {
      logger('create kmer file')
      identify_kmers(workDir, kmerSize, HLAfastaLoc, jellyFish)
    } else {
      logger('fishingStep', 'already done OR fishing turned off, skipping')
    }

    for (BAMfile in BAMfiles) {
      BAMid <- unlist(strsplit(BAMfile, split = '.bam'))[1]
      sample_dir <- paste(workDir, '/', BAMid, sep = '')
      if (!dir.exists(sample_dir)) {
        dir.create(sample_dir, recursive = TRUE)
      }

      ## extract HLA possible reads from BAM file
      logger('extract HLA possible reads from BAM file')

      ## chr 6 and contigs
      BAM_chr6_region_fn <- file.path(sample_dir,
        paste0(BAMid, '.chr6region.sam'))
      BAMfile_full <- file.path(BAMDir, BAMfile)
      ## copy header of original bam file
      samtoolsCMD <- paste0(samtools, ' view -H ', BAMfile_full, ' > ',
        BAM_chr6_region_fn)
      logger(samtoolsCMD)
      system(samtoolsCMD)

      if (!file.exists(BAM_chr6_region_fn) ||
        file.size(BAM_chr6_region_fn) <= 4) {
        warning(paste0('Did not correctly parse bam head for ', BAMfile_full))
      }
      header_file_size <- file.size(BAM_chr6_region_fn)

      # Read out the names of alternate loci
      # readLines(BAM_chr6_region_fn) %>%
      #   { grep('SN', ., value = T) } %>%
      #   { grep('Un', ., value = T) } %>%
      #   strsplit('\t') %>% .[[1]] %>% .[2] %>% { grepl('chr', .) }

      if (genomeAssembly == 'hg19') {
        hla_locs <- c('6:29909037-29913661', '6:31321649-31324964',
                      '6:31236526-31239869')
      } else if (genomeAssembly == 'grch38') {
        hla_locs <- c('6:29941260-29945884', '6:31353872-31357187',
                      '6:31268749-31272092')
      }

      ## Check if chromosome names are prefixed by the string 'chr'
      chr_prefix <- readLines(BAM_chr6_region_fn, n = 2)[2] %>%
        strsplit('\t') %>% .[[1]] %>% .[2] %>% { grepl('chr', .) }
      if (chr_prefix)
        hla_locs <- paste0('chr', hla_locs)

      for (i in seq_along(hla_locs)) {
        samtoolsCMD <- paste(samtools, ' view ', BAMfile_full, ' ',
          hla_locs[i], ' >> ', BAM_chr6_region_fn, sep = '')
        logger(samtoolsCMD)
        system(samtoolsCMD)
      }

      if (file.size(BAM_chr6_region_fn) <= header_file_size) {
        msg <- paste0('Did not filter reads correctly for ', BAMfile_full)
        howToWarn(msg)
        logger(msg)
        next
      }

      if (genomeAssembly == 'hg19') {
        alt_loci_names <- c('6_apd_hap1', '6_cox_hap2', '6_dbb_hap3',
          '6_mann_hap4', '6_mcf_hap5', '6_qbl_hap6', '6_ssto_hap7')
      } else if (genomeAssembly == 'grch38') {
        alt_loci_names <- c( 'GL000250.2', 'GL000251.2', 'GL000252.2',
          'GL000253.2', 'GL000254.2',  'GL000255.2', 'GL000256.2',  'KI270758.1',
          'NT_167244.2', 'NT_113891.3', 'NT_167245.2', 'NT_167246.2',
          'NT_167247.2', 'NT_167248.2', 'NT_167249.2', 'NT_187692.1')
      }
      if (chr_prefix && genomeAssembly == 'hg19')
        alt_loci_names <- paste0('chr', alt_loci_names)

      ## Strip away version numbers if they're present
      alt_loci_names <- alt_loci_names %>% { gsub('\\.\\d+', '', .) }

      ## Identify sequence names to which reads are mapped in the
      ## bam file
      sequence_names <- readLines(BAM_chr6_region_fn) %>%
        { .[grepl('@SQ', .)] } %>%
        { sapply(., function(x) strsplit(x, '\t')[[1]][2]) } %>%
        { gsub('SN:', '', .) } %>%
        { gsub('chrUn_', '', .) } %>%
        { gsub('v\\d$', '', .) } %>%
        setNames(NULL)

      alt_loci_names <- intersect(alt_loci_names, sequence_names)

      sam_file <- paste0(sample_dir, '/', BAMid, '.chr6region.sam')
      chr6.f1 <- paste0(sample_dir, '/', BAMid, '.chr6region.1.fastq')
      chr6.f2 <- paste0(sample_dir, '/', BAMid, '.chr6region.2.fastq')
      fished.f1 <- paste0(sample_dir, '/', 'fished.1.fastq')
      fished.f2 <- paste0(sample_dir,  '/', 'fished.2.fastq')
      hlaBAMfile <- paste0(sample_dir, '/', BAMid, '.chr6region.patient',
        '.reference.hlas.csorted.noduplicates.filtered.bam')

      if (length(alt_loci_names) > 0) {
        ## Add reads mapped to alternate loci to fished file
        for (chrom_name in alt_loci_names) {
          samtoolsCMD <- paste0(samtools, ' view ', BAMfile_full,
            ' ', chrom_name, ' >> ', sam_file)
          logger(samtoolsCMD)
          system(samtoolsCMD)
        }
      } else {
        logger('No alternate loci to fish in')
      }

      ## Turn into fastq -- this step has an error with unpaired mates, but
      ## seems to work ok just the same (VALIDATION_STRINGENCY = SILENT)
      logger('Turn (fished) reads into a fastq file')
      samToFastQ <- paste('java',
        ' -jar ', gatkDir, '/SamToFastq.jar',
        ' I=', sam_file,
        ' F=', chr6.f1,
        ' F2=', chr6.f2,
        ' VALIDATION_STRINGENCY=SILENT', sep = '')
      logger(samToFastQ)
      system(samToFastQ)

      # system(glue('wc -l {chr6.f1} {chr6.f2}'))

      if (fishingStep) {
        logger('Adding reads aligned to alternate loci')
        ## Generate fished.f1 & fished.f2 fastqs
        get.partially.matching.reads(workDir, sample_dir, BAMDir, BAMfile,
          kmerSize, gatkDir)
        logger('Combine chr6 reads with fished reads from alternate loci')
        ## Overwrite chr6.f1 and chr6.2 to include fished reads
        combine.fastqs(chr6.f1, chr6.f2, fished.f1, fished.f2)
      }

      logger(glue('{BAMfile}: aligning (fished) reads to all HLA alleles'))

      alignCMD <- paste0(novoDir, '/novoalign -d ', patient_hla_fasta_index_fn,
        ' -f ', chr6.f1, ' ', chr6.f2,
        ' -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW',
        ' -c ', novoThreads,
        ' 1> ', sample_dir, '/', BAMid,
        '.chr6region.patient.reference.hlas.sam',
        ' 2> ', sample_dir, '/', BAMid,
        '_BS_GL.chr6region.patient.reference.hlas.metrics')
      logger(alignCMD)
      system(alignCMD)

      ## Convert sam file to bam
      convertToBam <- paste(samtools, ' view -bS -o ',
        sample_dir, '/', BAMid, '.chr6region.patient.reference.hlas.bam', ' ',
        sample_dir, '/', BAMid, '.chr6region.patient.reference.hlas.sam',
        sep = '')
      logger(convertToBam)
      system(convertToBam)

      sortBAM <- paste('java -jar ', gatkDir, '/SortSam.jar',
        ' I=', sample_dir, '/', BAMid, '.chr6region.patient.reference.hlas.bam',
        ' O=', sample_dir, '/', BAMid,
        '.chr6region.patient.reference.hlas.csorted.bam',
        ' SORT_ORDER=coordinate', sep = '')
      logger(sortBAM)
      system(sortBAM)

      ## Remove duplicates
      removeDup <- paste(samtools, ' rmdup ',
        sample_dir, '/', BAMid,
        '.chr6region.patient.reference.hlas.csorted.bam',
        ' ',
        sample_dir, '/', BAMid,
        '.chr6region.patient.reference.hlas.csorted.noduplicates.bam', sep = '')
      logger(removeDup)
      system(removeDup)

      ## Only take reads that are in proper pair
      ## From samtools manual:
      ## -f INT
      ## Only output alignments with all bits set in INT present in the FLAG
      ## field. INT can be specified in hex by beginning with `0x' (i.e.
      ## /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/)
      ## [0].
      readPairs <- paste(samtools, ' view -f 2 -b -o ',
        hlaBAMfile,
        ' ',
        sample_dir, '/', BAMid,
        '.chr6region.patient.reference.hlas.csorted.noduplicates.bam',
        sep = '')
      logger(readPairs)
      system(readPairs)

      ## index the aligned bam
      indexBAM <- paste0(samtools, ' index ', hlaBAMfile)
      logger(indexBAM)
      system(indexBAM)

      for (allele in hlaAlleles) {
        logger(sprintf('get HLA specific SAM for allele: %s', allele))

        BAM_tmp_fn <- paste0(sample_dir, '/', BAMid, '.temp.', allele, '.bam')
        BAM_fn <- paste0(sample_dir, '/', BAMid, '.type.', allele, '.bam')
        BAM_filtered_fn <- paste0(sample_dir, '/', BAMid, '.type.', allele,
          '.filtered.bam')
        passed_reads_fn <- paste0(sample_dir, '/', BAMid, '.',
          allele, '.passed.reads.txt')

        getReads <- paste0(samtools, ' view -b -o ', ' ',
          BAM_tmp_fn, ' ', hlaBAMfile, ' ', allele)
        logger(getReads)
        system(getReads)

        if (is.null(BAM_tmp_fn) || is.na(BAM_tmp_fn) ||
          !file.exists(BAM_tmp_fn) || file.size(BAM_tmp_fn) <= 464) {
          msg <- sprintf('No reads mapped to %s', allele)
          howToWarn(msg)
          logger(msg)
          # system(sprintf('samtools view %s', BAM_filtered_fn))
          # system(sprintf('ls -ltr %s', BAM_filtered_fn))
        }

        samtoolsSort <- paste0(samtools, ' sort ', BAM_tmp_fn, ' -o ', BAM_fn)
        logger(samtoolsSort)
        system(samtoolsSort)

        samtoolsIndex <- paste0(samtools, ' index ', BAM_fn)
        logger(samtoolsIndex)
        system(samtoolsIndex)

        ## And filter out reads that have too many events -- has to be done here
        ## because some reads map to multiple alleles
        passed_reads <- count_events(BAM_fn, n = numMisMatch)
        if (is.null(passed_reads)) {
          msg <- sprintf('Could not get mapping reads for %s', allele)
          howToWarn(msg)
          logger(msg)
          # system(sprintf('cp %s %s', BAM_fn, BAM_filtered_fn))
        } else {
          write_tsv(passed_reads, passed_reads_fn, col.names = F)

          extractCMD <- paste0('java -jar ', gatkDir, '/FilterSamReads.jar',
            ' I=', BAM_fn,
            ' FILTER=includeReadList',
            ' READ_LIST_FILE=', passed_reads_fn,
            ' OUTPUT=', BAM_filtered_fn)
          logger(extractCMD)
          system(extractCMD)

          samtoolsIndex <- paste0(samtools, ' index ', BAM_filtered_fn)
          logger(samtoolsIndex)
          system(samtoolsIndex)
        }
      }
    }
  }
  ### }}}

  ## {{{ Get coverage for samples
  if (coverageStep) {
    for (sample in samples) {
      logger(sprintf('get coverage of HLA alleles for sample: %s', sample))
      sample_dir <- paste(workDir, '/', sample, sep = '')
      filteredBAMfiles <- grep('filtered.bam$',
        grep('type', list.files(sample_dir, recursive = T), value = TRUE),
        value = TRUE)
      if (length(filteredBAMfiles) == 0) {
        msg <- sprintf('No filtered bam files were made for %s, aborting LOHHLA',
          sample)
        logger(msg)
        error_msg <- glue('{error_msg};no_filtered_bam_present_{sample}')
        next
      }

      if (sample %in% tumorName) {
        type <- 'tumor'
      } else {
        type <- 'normal'
      }

      ## 2019-04-09 17:24 Verify not all filtered bam files are identical
      if (any(table(tools::md5sum(file.path(sample_dir, filteredBAMfiles))) > 1)) {
        logger('File sizes of filtered bam files: ',
          paste(file.size(file.path(sample_dir, filteredBAMfiles)),
            collapse = ', '))

        logger('MD5 sums of filtered bam files: ',
          paste(tools::md5sum(file.path(sample_dir, filteredBAMfiles)),
            collapse = ', '))
      }
      # file.size(filteredBAMfiles)
      # file.exists(filteredBAMfiles)

      ## Get pileup files for each bam
      for (BAMfile in filteredBAMfiles) {
        hlaAllele <- grep(pattern = 'hla',
          x = unlist(strsplit(BAMfile, split = '\\.')), value = TRUE)
        ## Make hlaAllele inference robust to file names that have the words
        #$ hla_loci in them
        if (length(hlaAllele) > 1) {
          hlaAllele <- grep('hla_(?!loci)', hlaAllele, value = T, perl = T)
        }
        mpileupFile <- paste0(workDir, '/', sample, '.',
          hlaAllele, '.', type, '.mpileup')

        cmd <- paste(samtools, ' mpileup ', sample_dir, '/', BAMfile, ' -f ',
                     HLAfastaLoc, ' > ', mpileupFile, sep = '')

        if (file.exists(mpileupFile)) {
          logger(paste(mpileupFile, 'already exists, skipping.'))
        } else {
          logger(cmd)
          system(cmd)
        }
        logger(system(glue('head {mpileupFile} -n 1'), intern = T))
      }
    }
    test_error_presence <- sapply(samples, function(sample) {
      grepl(glue('no_filtered_bam_present_{sample}'), error_msg)
    }) %>% any(na.rm = T)
    if (test_error_presence) {
      coverageStep <- F
    }
  }
  ### }}}

  ## {{{ Extract number of unique reads sequenced in tumor and normal
  if (any(unique(runWithNormal)) && 
      (is.null(opt$normalAlignedReads) || is.null(opt$tumorAlignedReads))) {
    if (!override) {
      sample_uniq_mapped_reads <-
        getUniqMapReads(workDir = workDir, BAMDir = BAMDir,
          tumorBAMfile = tumorBAMfile, normalBAMfile = normalBAMfile,
          override = FALSE)
    } else {
      sample_uniq_mapped_reads <- getUniqMapReads(workDir = workDir,
        BAMDir = BAMDir, override = TRUE, overrideDir = overrideDir)
    }

    ## In case of multiple normal samples, reduce to the best normal sample
    ## based on maximal read coverage
    normal_sample <- unlist(sample_uniq_mapped_reads)[normalName]
    if (length(normal_sample) > 1) {
      best_normal_sample <- names(normal_sample[which.max(normal_sample)])
      normalName <- best_normal_sample
      samples <- samples[samples %in% c(tumorName, normalName)]
    }

    ## In case of multiple tumor samples, reduce to the best tumor sample
    tumor_sample <- unlist(sample_uniq_mapped_reads)[tumorName]
    if (length(tumor_sample) > 1) {
      best_tumor_sample <- names(tumor_sample[which.max(tumor_sample)])
      tumorName <- best_tumor_sample
      samples <- samples[samples %in% c(tumorName, normalName)]
    }

    GermLineUniqMappedReads <- unlist(sample_uniq_mapped_reads[normalName])
  } else {
    GermLineUniqMappedReads <- opt$normalAlignedReads
  }
  ### }}}

  ## {{{ Compare coverage between alleles
  ## Load the winners.  Next, we can look at each mpileupFile, and assess
  ## whether we see differences in coverage between the two.  Look at a
  ## sample of interest.
  ## Don't process normal samples
  for (sample in tumorName) {
    if (runWithNormal) {
      if (is.null(opt$tumorAlignedReads)) {
        UniqMappedReads <- sample_uniq_mapped_reads[[sample]]
      } else {
        UniqMappedReads <- opt$tumorAlignedReads
      }
      MultFactor <- as.numeric(GermLineUniqMappedReads) /
        as.numeric(UniqMappedReads)
    } else {
      MultFactor <- 1
    }

    combinedTable <- NULL
    hlas <- c('hla_a', 'hla_b', 'hla_c')
    HLAoutPut_l <- purrr::map(hlas, function(HLA_gene) {
      HLA_As <- grep(HLA_gene, hlaAlleles, value = TRUE)
      ## Change this to a test of whether coverage files are present
      if (!coverageStep) {
        return(list(message = glue('did_not_perform_coverage_{error_msg}'),
            HLA_A_type1 = repl_NA(HLA_As[1]),
            HLA_A_type2 = repl_NA(HLA_As[2])))
      } else {
        if (F && HLA_gene == 'hla_c') browser()
        # browser()
        logger(sprintf('analyzing coverage differences in sample: %s, hla: %s',
            sample, HLA_gene))

        region_HLA_plot_data <- paste(figureDir, '/', sample, '.', HLA_gene,
          '.tmp.data.plots.RData', sep = '')
        ## Remove tmp plotting files so if this step fails, can't plot
        ## incorrectly
        if (file.exists(region_HLA_plot_data)) {
          cmd <- sprintf('rm %s', region_HLA_plot_data)
          logger(cmd)
          system(cmd)
        }

        if (length(HLA_As) == 1) {
          return(list(message = 'homozygous_alleles_not_implemented',
              HLA_A_type1 = repl_NA(HLA_As[1]),
              HLA_A_type2 = repl_NA(HLA_As[2])))
        } else if (length(HLA_As) == 0) {
          return(list(message = 'no_recognized_alleles'))
        }
        HLA_A_type1 <- HLA_As[1]
        HLA_A_type2 <- HLA_As[2]

        ## The reference sequence for the patient's two HLA alleles
        HLA_type1Fasta <- hlaFasta[[HLA_A_type1]]
        HLA_type2Fasta <- hlaFasta[[HLA_A_type2]]

        ## Perform local pairwise alignment
        seqs <- c(PasteVector(toupper(HLA_type1Fasta), sep = ""),
                  PasteVector(toupper(HLA_type2Fasta), sep = ""))
        sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1,
          baseOnly = TRUE)

        missMatchPositions <- pairwiseAlignment(seqs[1], seqs[2],
          substitutionMatrix = sigma,
          gapOpening = -2, gapExtension = -4,
          scoreOnly = FALSE, type = 'local') %>%
          getMisMatchPositionsPairwiseAlignment(returnlist = TRUE)

        if (length(missMatchPositions$diffSeq1) == 0) {
          msg <- 'HLA alleles are identical, although they have different names'
          logger(msg)
          return(list(message = 'alleles_identical',
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }
        if (length(missMatchPositions$diffSeq1) < 5) {
          msg <- glue::glue('HLA alleles are very similar (fewer than 5 \\
            mismatch positions)! Keep that in mind when considering results.')
          logger(msg)
          warning(msg)
        }

        HLA_A_type1normalLoc <- grep(pattern = HLA_A_type1,
          x = list.files(workDir, pattern = 'normal\\.mpileup$',
            full.names = TRUE), value = TRUE)

        HLA_A_type2normalLoc <- grep(pattern = HLA_A_type2,
          x = list.files(workDir, pattern = "normal\\.mpileup$",
            full.names = TRUE), value = TRUE)

        if (length(HLA_A_type1normalLoc) == 0 ||
            length(HLA_A_type2normalLoc) == 0 ||
            file.size(HLA_A_type1normalLoc) == 0 ||
            file.size(HLA_A_type2normalLoc) == 0) {
          msg <- sprintf('No coverage in normal sample, aborting %s', HLA_gene)
          howToWarn(msg)
          logger(msg)
          return(list(message = 'no_pileup_for_normal',
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }

        if (any(runWithNormal)) {
          HLA_A_type1normal <- tryCatch(read.table(HLA_A_type1normalLoc,
            sep = '\t', stringsAsFactors = FALSE, quote = '', fill = TRUE), 
          error = function(e) { logger(
            glue('Could not access {HLA_A_type1normalLoc}')) }) 
          if (is.null(HLA_A_type1normal)) {
            return(list(message = glue('no_pileup_for_normal_{HLA_gene}_1'),
                HLA_A_type1 = HLA_A_type1,
                HLA_A_type2 = HLA_A_type2))
          }

          HLA_A_type2normal <- tryCatch(read.table(HLA_A_type2normalLoc,
            sep = '\t', stringsAsFactors = FALSE, quote = '', fill = TRUE), 
          error = function(e) { logger(
            glue('Could not access {HLA_A_type2normalLoc}')) }) 
          if (is.null(HLA_A_type2normal)) {
            return(list(message = glue('no_pileup_for_normal_{HLA_gene}_2'),
                HLA_A_type1 = HLA_A_type1,
                HLA_A_type2 = HLA_A_type2))
          }
        } else {
          ## {{{ Type 1
          HLA_A_type1normal <- data.frame(
            cbind(HLA_A_type1, 1:length(HLA_type1Fasta),
              toupper(HLA_type1Fasta), minCoverageFilter + 1),
            stringsAsFactors = FALSE)
          colnames(HLA_A_type1normal) <-
            paste('V', 1:ncol(HLA_A_type1normal), sep = '')
          HLA_A_type1normal$V4 <- as.numeric(HLA_A_type1normal$V4)
          ## }}} Type 1
          ## {{{ Type 2
          HLA_A_type2normal <- data.frame(
            cbind(HLA_A_type2, 1:length(HLA_type2Fasta),
              toupper(HLA_type2Fasta), minCoverageFilter + 1),
            stringsAsFactors = FALSE)
          colnames(HLA_A_type2normal) <-
            paste('V', 1:ncol(HLA_A_type2normal), sep = '')
          HLA_A_type2normal$V4 <- as.numeric(HLA_A_type2normal$V4)
          ## }}} Type 2
        }

        ## 2019-04-04 11:35 M.S. Rownames might not be unique, using the
        ## rownames attribute is therefore unsafe
        ## Type 1 {{{
        HLA_A_type1normal$mm_position <- as.character(HLA_A_type1normal$V2)
        tumpile <- paste0(workDir, '/', sample, '.', 
          HLA_A_type1, '.', 'tumor.mpileup')
        HLA_A_type1tumor <- tryCatch(
          read.table(tumpile, sep = '\t', stringsAsFactors = FALSE, 
            quote = '', fill = TRUE, col.names = paste0('V', c(1:6))), 
          error = function(e) { print(e); NULL }) 
        if (is.null(HLA_A_type1tumor)) {
          return(list(message = glue('no_pileup_for_tumor_{HLA_gene}_1'),
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }

        HLA_A_type1tumor$mm_position <- HLA_A_type1tumor$V2
        ## Apply minimum coverage thresholds (we only apply this to the normal
        ## for now)
        HLA_A_type1normal <-
          HLA_A_type1normal[HLA_A_type1normal$V4 > minCoverageFilter,
          , drop = FALSE]
        ## Select loci for which coverage info is available for both tumor and
        ## matched normal
        tmp <- intersect(HLA_A_type1tumor$mm_position,
          HLA_A_type1normal$mm_position)
        HLA_A_type1normal <-
          HLA_A_type1normal[HLA_A_type1normal$mm_position %in% tmp, , drop = FALSE] %>% unique
        # HLA_A_type1normal[duplicated(HLA_A_type1normal$mm_position), ]
        # HLA_A_type1normal[HLA_A_type1normal$mm_position == '220', ]
        # HLA_A_type1normal[HLA_A_type1normal$mm_position == '247', ]
        # HLA_A_type1tumor[HLA_A_type1tumor$mm_position == '247', ]

        HLA_A_type1tumor <-
          HLA_A_type1tumor[HLA_A_type1tumor$mm_position %in% tmp, , drop = FALSE] %>% unique

        HLA_A_type1normalCov <- HLA_A_type1normal$V4
        names(HLA_A_type1normalCov) <- HLA_A_type1normal$V2
        HLA_A_type1normalCov <-
          HLA_A_type1normalCov[HLA_A_type1tumor$mm_position]

        HLA_A_type1tumorCov <- rep(0, length(HLA_A_type1normalCov))
        names(HLA_A_type1tumorCov) <- names(HLA_A_type1normalCov)
        HLA_A_type1tumorCov[HLA_A_type1tumor$mm_position] <- HLA_A_type1tumor$V4
        ## }}} Type 1

        ## Type 2 {{{
        HLA_A_type2normal$mm_position <- as.character(HLA_A_type2normal$V2)
        tumpile <- paste0(workDir, '/', sample, '.', HLA_A_type2, 
          '.', 'tumor.mpileup')
        HLA_A_type2tumor <- tryCatch(
          read.table(tumpile, sep = '\t', stringsAsFactors = FALSE, 
            quote = '', fill = TRUE, col.names = paste0('V', c(1:6))), 
          error = function(e) { print(e); NULL }) 
        if (is.null(HLA_A_type2tumor)) {
          return(list(message = glue('no_pileup_for_tumor_{HLA_gene}_2'),
              HLA_A_type2 = HLA_A_type2,
              HLA_A_type2 = HLA_A_type2))
        }

        HLA_A_type2tumor$mm_position <- HLA_A_type2tumor$V2
        # if (F && HLA_gene == 'hla_b') browser() 
        ## Apply minimum coverage thresholds (we only apply this to the normal
        ## for now)
        HLA_A_type2normal <-
          HLA_A_type2normal[HLA_A_type2normal$V4 > minCoverageFilter,
          , drop = FALSE]
        ## Select loci for which coverage info is available for both tumor and
        ## matched normal
        tmp <- intersect(HLA_A_type2tumor$mm_position,
          HLA_A_type2normal$mm_position)
        HLA_A_type2normal <-
          HLA_A_type2normal[HLA_A_type2normal$mm_position %in% tmp, , drop = FALSE] %>% unique

        HLA_A_type2tumor <-
          HLA_A_type2tumor[HLA_A_type2tumor$mm_position %in% tmp, , drop = FALSE] %>% unique

        HLA_A_type2normalCov <- HLA_A_type2normal$V4
        names(HLA_A_type2normalCov) <- HLA_A_type2normal$V2
        HLA_A_type2normalCov <-
          HLA_A_type2normalCov[HLA_A_type2tumor$mm_position]

        HLA_A_type2tumorCov <- rep(0, length(HLA_A_type2normalCov))
        names(HLA_A_type2tumorCov) <- names(HLA_A_type2normalCov)
        HLA_A_type2tumorCov[HLA_A_type2tumor$mm_position] <- HLA_A_type2tumor$V4
        ## }}} Type 2

        ## catch issues with HLA coverage
        HLA_type1_ok <- intersect(names(HLA_A_type1normalCov),
          missMatchPositions$diffSeq1) %>% length
        HLA_type2_ok <- intersect(names(HLA_A_type2normalCov),
          missMatchPositions$diffSeq2) %>% length

        if (nrow(HLA_A_type1normal) == 0 || nrow(HLA_A_type2normal) == 0 ||
          HLA_type1_ok == 0 || HLA_type2_ok == 0) {
          msg <- glue::glue('No position has greater than minimum \\
            coverage filter for {HLA_gene}')
          logger(msg)
          return(list(message = 'no_position_sufficiently_covered_in_normal',
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }

        if (HLA_type1_ok / HLA_type2_ok < 0.05) {
          msg <- paste('Check that the HLA type is correct for ',
            HLA_A_type1, '!', sep = '')
          logger(msg)
          howToWarn(msg)
          return(list(
              message = glue('mismatch_pos_density_imbalance_in_normal_allele1'),
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }
        if (HLA_type2_ok / HLA_type1_ok < 0.05) {
          msg <- paste('Check that the HLA type is correct for ',
            HLA_A_type2, '!', sep = '')
          logger(msg)
          howToWarn(msg)
          return(list(
              message = glue('mismatch_pos_density_imbalance_in_normal_allele2'),
              HLA_A_type1 = HLA_A_type1,
              HLA_A_type2 = HLA_A_type2))
        }

        if (extractNONmismatchReads == T) {
          ## Next, let's extract the reads that do not cover any mismatches,
          ## and then mpileup them.
          ## First, make bed-files for mismatch positions
          ## {{{ Type 1
          mismatchPosSeq1 <-
            HLA_A_type1tumor[HLA_A_type1tumor$V2 %in%
            missMatchPositions$diffSeq1, , drop = FALSE]
          missMatchBed1 <- mismatchPosSeq1[, c(1, 2, 2)]
          ## 2019-05-01 11:00 M.S. Changed the bed file from three to two
          ## columns
          missMatchBed1$V2 <- as.numeric(missMatchBed1$V2) - 1
          missMatchBed1$V2.1 <- as.numeric(missMatchBed1$V2.1) + 1
          type1_bed <- paste(workDir, '/', sample, '.', HLA_A_type1, '.bed',
            sep = '')
          write.table(missMatchBed1,
            file = type1_bed,
            quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)
          ## }}} Type 1
          ## {{{ Type 2
          mismatchPosSeq2 <-
            HLA_A_type2tumor[HLA_A_type2tumor$V2 %in%
            missMatchPositions$diffSeq2, , drop = FALSE]
          missMatchBed2 <- mismatchPosSeq2[, c(1, 2, 2)]
          missMatchBed2$V2 <- as.numeric(missMatchBed2$V2) - 1
          missMatchBed2$V2.1 <- as.numeric(missMatchBed2$V2.1) + 1
          type2_bed <- paste(workDir, '/', sample, '.', HLA_A_type2, '.bed',
            sep = '')
          write.table(missMatchBed2,
            file = type2_bed,
            quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)
          ## }}} Type 2

          ## {{{ Type 1
          ## Select regions that are not covered in bed file
          Type1TumorCmd <- paste(bedtools, ' intersect -v ',
            ' -a ', workDir, '/', sample, '/', sample, '.type.', HLA_A_type1,
            '.filtered.bam',
            ' -b ', workDir, '/', sample, '.', HLA_A_type1, '.bed',
            ' > ',
            workDir, '/', sample, '.', HLA_A_type1, '.tumor.NoMissMatch.bam',
            sep = '')
          logger(Type1TumorCmd)
          system(Type1TumorCmd)

          if (runWithNormal) {
            Type1NormalCmd <- paste(bedtools, ' intersect -v ',
              '-a ', workDir, '/', normalName, '/', normalName,
              '.type.', HLA_A_type1, '.filtered.bam ',
              '-b ', workDir, '/', sample, '.', HLA_A_type1, '.bed',
              ' > ',
              workDir, '/', sample, '.', HLA_A_type1, '.normal.NoMissMatch.bam',
              sep = '')
            system(Type1NormalCmd)
          }
          ## }}} Type 1
          ## {{{ Type 2
          Type2TumorCmd <- paste(bedtools, " intersect -v -a ", workDir, '/',
            sample, "/", sample, ".type.", HLA_A_type2, ".filtered.bam",
            " -b ", workDir, '/', sample, ".", HLA_A_type2, ".bed",
            " > ",
            workDir, '/', sample, ".", HLA_A_type2, ".tumor.NoMissMatch.bam",
            sep = "")
          logger(Type2TumorCmd)
          system(Type2TumorCmd)

          if (any(runWithNormal)) {
            Type2NormalCmd <- paste(bedtools, ' intersect -v -a ',
              workDir, '/', normalName, '/', normalName, '.type.', HLA_A_type2,
              '.filtered.bam -b ',
              workDir, '/', sample, '.', HLA_A_type2, '.bed', ' > ',
              workDir, '/', sample, '.', HLA_A_type2, '.normal.NoMissMatch.bam',
              sep = '')
            system(Type2NormalCmd)
          }
          ## }}} Type 2

          ## next, let's do the mpileup step.
          ## {{{ Type 1
          nomismatch_bam <- paste0(workDir, '/', sample, '.',
            HLA_A_type1, '.tumor.NoMissMatch.bam')
          nomismatch_pu <- paste0(workDir, '/', sample, '.',
            HLA_A_type1, '.tumor.NoMissMatch.pileup')
          MpilupType1TumorCmd <- paste(samtools, ' mpileup ',
            nomismatch_bam, ' -f ', HLAfastaLoc, ' > ', nomismatch_pu, sep = '')
          logger(MpilupType1TumorCmd)
          system(MpilupType1TumorCmd)

          if (any(runWithNormal)) {
            nomismatch_bam <- paste0(workDir, '/', sample, '.',
              HLA_A_type1, '.normal.NoMissMatch.bam')
            nomismatch_pu <- paste0(workDir, '/', sample, '.',
              HLA_A_type1, '.normal.NoMissMatch.pileup')

            MpilupType1NormalCmd <- paste(samtools, ' mpileup ',
              nomismatch_bam, ' -f ', HLAfastaLoc, ' > ', nomismatch_pu,
              sep = '')
            logger(MpilupType1NormalCmd)
            system(MpilupType1NormalCmd)
          }
          ## }}} Type 1
          ## {{{ Type 2
          nomismatch_bam <- paste0(workDir, '/', sample, '.',
            HLA_A_type2, '.tumor.NoMissMatch.bam')
          nomismatch_pu <- paste0(workDir, '/', sample, '.',
            HLA_A_type2, '.tumor.NoMissMatch.pileup')
          MpilupType2TumorCmd <- paste(samtools, ' mpileup ',
            nomismatch_bam, ' -f ', HLAfastaLoc, ' > ', nomismatch_pu, sep = '')
          logger(MpilupType2TumorCmd)
          system(MpilupType2TumorCmd)

          if (any(runWithNormal)) {
            nomismatch_bam <- paste0(workDir, '/', sample, '.',
              HLA_A_type2, '.normal.NoMissMatch.bam')
            nomismatch_pu <- paste0(workDir, '/', sample, '.',
              HLA_A_type2, '.normal.NoMissMatch.pileup')

            MpilupType2NormalCmd <- paste(samtools, ' mpileup ',
              nomismatch_bam, ' -f ', HLAfastaLoc, ' > ', nomismatch_pu,
              sep = '')
            logger(MpilupType2NormalCmd)
            system(MpilupType2NormalCmd)
          }
          ## }}} Type 2

          ## Get the coverage for the sites
          ## {{{ Type 1
          nomismatch_pu <- paste(workDir, '/', sample, ".",
            HLA_A_type1, ".tumor.NoMissMatch.pileup", sep = "")
          HLA_type1tumor_nomissmatch <- tryCatch(read.table(nomismatch_pu,
              stringsAsFactors = FALSE, fill = TRUE, quote = "", sep = '\t'),
            error = function(e) NULL)
          ## }}} Type 1
          ## {{{ Type 2
          nomismatch_pu <- paste(workDir, '/', sample, ".",
            HLA_A_type2, ".tumor.NoMissMatch.pileup", sep = "")
          HLA_type2tumor_nomissmatch <- tryCatch(read.table(nomismatch_pu,
              stringsAsFactors = FALSE, fill = TRUE, quote = "", sep = '\t'),
            error = function(e) NULL)
          ## }}} Type 2

          if (runWithNormal) {
            ## {{{ Type 1
            nomismatch_pu <- paste(workDir, '/', sample, ".",
              HLA_A_type1, ".normal.NoMissMatch.pileup", sep = "")
            HLA_type1normal_nomissmatch <- tryCatch(read.table(nomismatch_pu,
                stringsAsFactors = FALSE, fill = TRUE, quote = "", sep = '\t'),
              error = function(e) NULL)
            ## }}} Type 1
            ## {{{ Type 2
            nomismatch_pu <- paste(workDir, '/', sample, ".",
              HLA_A_type2, ".normal.NoMissMatch.pileup", sep = "")
            HLA_type2normal_nomissmatch <- tryCatch(read.table(nomismatch_pu,
                stringsAsFactors = FALSE, fill = TRUE, quote = "", sep = '\t'),
              error = function(e) NULL)
            ## }}} Type 2
          } else {
            ## Make a dummy one if there isn't a normal
            ## {{{ Type 1
            HLA_type1normal_nomissmatch <-
              data.frame(cbind(HLA_A_type1, 1:length(HLA_type1Fasta),
                  toupper(HLA_type1Fasta), minCoverageFilter + 1),
                stringsAsFactors = FALSE)
            HLA_type1normal_nomissmatch$V4 <-
              as.numeric(HLA_type1normal_nomissmatch$V4)
            ## }}} Type 1
            ## {{{ Type 2
            HLA_type2normal_nomissmatch <-
              data.frame(cbind(HLA_A_type2, 2:length(HLA_type2Fasta),
                  toupper(HLA_type2Fasta), minCoverageFilter + 2),
                stringsAsFactors = FALSE)
            HLA_type2normal_nomissmatch$V4 <-
              as.numeric(HLA_type2normal_nomissmatch$V4)
            ## }}} Type 2
          }

          ## error if non-mismatch file is empty
          ## can probably work around missing non-mismatches later...
          if (is.null(HLA_type1tumor_nomissmatch) ||
              is.null(HLA_type2tumor_nomissmatch) ||
              is.null(HLA_type1normal_nomissmatch) ||
              is.null(HLA_type2normal_nomissmatch)) {
            msg <- paste('No non-mismatch positions for ', HLA_gene, sep = '')
            logger(msg)
            # return(list(
            #     message = glue('no_non_mismatch_positions_{HLA_gene}'),
            #     HLA_A_type1 = HLA_A_type1,
            #     HLA_A_type2 = HLA_A_type2))
          } else {
            ## Apply minimum coverage thresholds (we only apply this to the
            ## normal for now)
            ## {{{ Type 1
            HLA_type1tumor_nomissmatch$rownames <-
              HLA_type1tumor_nomissmatch$V2
            HLA_type1normal_nomissmatch$rownames <-
              HLA_type1normal_nomissmatch$V2

            HLA_type1normal_nomissmatch %<>%
              { .[.$V4 > minCoverageFilter, , drop = FALSE] }

            tmp <- intersect(
              HLA_type1tumor_nomissmatch$rownames,
              HLA_type1normal_nomissmatch$rownames)
            HLA_type1tumor_nomissmatch %<>% dplyr::filter(rownames %in% tmp)
            HLA_type1normal_nomissmatchCov <- HLA_type1normal_nomissmatch$V4
            names(HLA_type1normal_nomissmatchCov) <- HLA_type1normal_nomissmatch$V2

            HLA_type1tumor_nomissmatchCov <-
              rep(0, length(HLA_type1normal_nomissmatchCov))
            names(HLA_type1tumor_nomissmatchCov) <-
              names(HLA_type1normal_nomissmatchCov)
            HLA_type1tumor_nomissmatchCov[HLA_type1tumor_nomissmatch$rownames] <-
              HLA_type1tumor_nomissmatch$V4
            ## }}} Type 1
            ## {{{ Type 2
            HLA_type2tumor_nomissmatch$rownames <-
              HLA_type2tumor_nomissmatch$V2
            HLA_type2normal_nomissmatch$rownames <-
              HLA_type2normal_nomissmatch$V2

            HLA_type2normal_nomissmatch <-
              HLA_type2normal_nomissmatch[HLA_type2normal_nomissmatch$V4>minCoverageFilter,
              , drop = FALSE]

            tmp <- intersect(
              HLA_type2tumor_nomissmatch$rownames,
              HLA_type2normal_nomissmatch$rownames)
            HLA_type2tumor_nomissmatch %<>% dplyr::filter(rownames %in% tmp)
            HLA_type2normal_nomissmatchCov <- HLA_type2normal_nomissmatch$V4
            names(HLA_type2normal_nomissmatchCov) <- HLA_type2normal_nomissmatch$V2

            HLA_type2tumor_nomissmatchCov <-
              rep(0, length(HLA_type2normal_nomissmatchCov))
            names(HLA_type2tumor_nomissmatchCov) <-
              names(HLA_type2normal_nomissmatchCov)
            HLA_type2tumor_nomissmatchCov[rownames(HLA_type2tumor_nomissmatch)] <-
              HLA_type2tumor_nomissmatch$V4
            ## }}} Type 2
          }
        }

        if (extractUniqueReads == TRUE) {
          ## Use the same mismatch positions, already have the bed to get reads
          ## that do overlap a mismatch
          ## {{{ type 1
          Type1TumorCmd <- paste(bedtools, ' intersect -loj -bed -b ',
            workDir, '/', sample, '/', sample, '.type.', HLA_A_type1, '.filtered.bam',
            ' -a ',
            workDir, '/', sample, '.', HLA_A_type1, '.bed', ' > ',
            workDir, '/', sample, '.', HLA_A_type1, '.tumor.mismatch.reads.bed',
            sep = '')
          logger(Type1TumorCmd)
          system(Type1TumorCmd)

          if (runWithNormal) {
            Type1NormalCmd <- paste(bedtools, ' intersect -loj -bed -b ',
              workDir, '/', normalName, '/', normalName, '.type.', HLA_A_type1,
              '.filtered.bam',
              ' -a ',
              workDir, '/', sample, '.', HLA_A_type1, '.bed', ' > ',
              workDir, '/', sample, '.', HLA_A_type1, '.normal.mismatch.reads.bed',
              sep = '')
            logger(Type1NormalCmd)
            system(Type1NormalCmd)
          }
          ## }}} type 1
          ## {{{ type 2
          Type2TumorCmd <- paste(bedtools, ' intersect -loj -bed -b ',
            workDir, '/', sample, '/', sample, '.type.', HLA_A_type2, '.filtered.bam',
            ' -a ',
            workDir, '/', sample, '.', HLA_A_type2, '.bed', ' > ',
            workDir, '/', sample, '.', HLA_A_type2, '.tumor.mismatch.reads.bed',
            sep = '')
          logger(Type2TumorCmd)
          system(Type2TumorCmd)

          if (runWithNormal) {
            Type2NormalCmd <- paste(bedtools, ' intersect -loj -bed -b ',
              workDir, '/', normalName, '/', normalName, '.type.', HLA_A_type2,
              '.filtered.bam',
              ' -a ',
              workDir, '/', sample, '.', HLA_A_type2, '.bed', ' > ',
              workDir, '/', sample, '.', HLA_A_type2, '.normal.mismatch.reads.bed',
              sep = '')
            logger(Type2NormalCmd)
            system(Type2NormalCmd)
          }
          ## }}} type 2
          ## Only take unique reads {{{
          hla_beds <- grep(pattern = sample,
            x = list.files(workDir, pattern = 'mismatch.reads.bed',
              full.names = TRUE), value = TRUE) %>%
            { grep(pattern = HLA_gene, x = ., value = TRUE) }
          for (i in hla_beds) {
            if (file.size(i) > 0) {
              x <- read.csv(i, sep = '\t', as.is = TRUE, header = FALSE)
              x <- x[!is.na(x$V2), ]
              x <- x[!duplicated(x$V8), ]
            } else {
              x <- t(c(paste('V', seq(1, 10), sep = '')))
            }
            write.table(x, file = gsub(pattern = 'mismatch.reads.bed',
                replacement = 'mismatch.unique.reads.bed', x = i),
              sep = '\t', quote = FALSE, row.names = FALSE)
          }
          ## }}}
          ## {{{ Type 1
          HLA_A_type1normalCov_mismatch <- HLA_A_type1normalCov %>%
            .[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1]
          HLA_A_type1normalCov_mismatch_unique <- HLA_A_type1normalCov_mismatch

          if (runWithNormal) {
            HLA_A_type1normal_unique <- read.table(paste(
                workDir, '/', sample, '.', HLA_A_type1,
                '.normal.mismatch.unique.reads.bed', sep = ''),
              sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
            HLA_A_type1normalCov_mismatch_unique <-
              table(HLA_A_type1normal_unique$V3) %>%
              .[names(HLA_A_type1normalCov_mismatch_unique)] %>%
              as.numeric %>%
              setNames(names(HLA_A_type1normalCov_mismatch)) %>%
              repl_NA(0)
            if (length(HLA_A_type1normalCov_mismatch_unique) != 0) {
              HLA_A_type1normalCov_mismatch_unique <-
                HLA_A_type1normalCov_mismatch_unique + 1
            }
          }

          HLA_A_type1tumorCov_mismatch <- HLA_A_type1tumorCov %>%
            .[names(HLA_A_type1tumorCov) %in% missMatchPositions$diffSeq1]
          HLA_A_type1tumor_unique <- read.table(
            paste(workDir, '/', sample, '.', HLA_A_type1,
              '.tumor.mismatch.unique.reads.bed', sep = ''),
            sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
          HLA_A_type1tumorCov_mismatch_unique <-
            table(HLA_A_type1tumor_unique$V3) %>%
            .[names(HLA_A_type1tumorCov_mismatch)] %>%
            as.numeric %>%
            setNames(names(HLA_A_type1tumorCov_mismatch)) %>%
            repl_NA(0)
          if (length(HLA_A_type1tumorCov_mismatch_unique) != 0) {
            HLA_A_type1tumorCov_mismatch_unique <-
              HLA_A_type1tumorCov_mismatch_unique + 1
          }
          ## }}} Type 1
          ## {{{ Type 2
          HLA_A_type2normalCov_mismatch <- HLA_A_type2normalCov %>%
            .[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2]
          HLA_A_type2normalCov_mismatch_unique <- HLA_A_type2normalCov_mismatch

          if (runWithNormal) {
            HLA_A_type2normal_unique <- read.table(paste(
                workDir, '/', sample, '.', HLA_A_type2,
                '.normal.mismatch.unique.reads.bed', sep = ''),
              sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
            HLA_A_type2normalCov_mismatch_unique <-
              table(HLA_A_type2normal_unique$V3) %>%
              .[names(HLA_A_type2normalCov_mismatch)] %>%
              as.numeric %>%
              setNames(names(HLA_A_type2normalCov_mismatch)) %>%
              repl_NA(0)
            if (length(HLA_A_type2normalCov_mismatch_unique) != 0) {
              HLA_A_type2normalCov_mismatch_unique <-
                HLA_A_type2normalCov_mismatch_unique + 1
            }
          }

          HLA_A_type2tumorCov_mismatch <- HLA_A_type2tumorCov %>%
            .[names(HLA_A_type2tumorCov) %in% missMatchPositions$diffSeq2]
          HLA_A_type2tumor_unique <- read.table(
            paste(workDir, '/', sample, '.', HLA_A_type2,
              '.tumor.mismatch.unique.reads.bed', sep = ''),
            sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
          HLA_A_type2tumorCov_mismatch_unique <-
            table(HLA_A_type2tumor_unique$V3) %>%
            .[names(HLA_A_type2tumorCov_mismatch)] %>%
            as.numeric %>%
            setNames(names(HLA_A_type2tumorCov_mismatch)) %>%
            repl_NA(0)
          if (length(HLA_A_type2tumorCov_mismatch_unique) != 0) {
            HLA_A_type2tumorCov_mismatch_unique <-
              HLA_A_type2tumorCov_mismatch_unique + 1
          }
          ## }}} Type 2
        }

        if (performIntegerCopyNum) {
          performIntegerCopyNumTmp <- performIntegerCopyNum
          copyNumSolutions <- read.table(CopyNumLoc, header = TRUE,
            stringsAsFactors = FALSE)
          colnames(copyNumSolutions) <- tolower(colnames(copyNumSolutions))

          if (!'tumorpurity' %in% colnames(copyNumSolutions) ||
              !'tumorploidy' %in% colnames(copyNumSolutions)) {
            msg <- glue::glue('column names tumorPloidy and tumorPurity are\\
              needed within {CopyNumLoc} if you wish to perform integer copy\\
              number')
            howToWarn(msg)
            logger(msg)
            performIntegerCopyNumTmp <- FALSE
          }

          if (!sample %in% rownames(copyNumSolutions)) {
            msg <- glue::glue('row names of {CopyNumLoc} must match sample\\
              names if you wish to perform integer copy number. \\
              Skipping sample: {sample}')
            howToWarn(msg)
            logger(msg)
            performIntegerCopyNumTmp <- FALSE
          }

          if (performIntegerCopyNumTmp) {
            tumorPloidy <- copyNumSolutions[sample, 'tumorploidy']
            tumorPurity <- copyNumSolutions[sample, 'tumorpurity']
          } else {
            tumorPloidy <- NA
            tumorPurity <- NA
          }
        }

        ## Infer copy number using combined BAF and logR
        misMatchCoveredInBoth <- cbind(
          ifelse(missMatchPositions$diffSeq1 %in% names(HLA_A_type1normalCov),
            1, 0),
          ifelse(missMatchPositions$diffSeq2 %in% names(HLA_A_type2normalCov),
            1, 0)) %>%
          { rowSums(.) == 2 }

        missMatchseq1 <- missMatchPositions$diffSeq1[misMatchCoveredInBoth]
        missMatchseq2 <- missMatchPositions$diffSeq2[misMatchCoveredInBoth]

        ## Let's get bins to collect for the coverage estimates
        startChar <- min(as.numeric(c(names(HLA_A_type1tumorCov),
          names(HLA_A_type2tumorCov))), na.rm = T)
        endChar <- max(as.numeric(c(names(HLA_A_type1tumorCov),
          names(HLA_A_type2tumorCov))), na.rm = T)
        seqToConsider <- seq(startChar, endChar, by = binSize)
        seqToConsider <- c(seqToConsider[-length(seqToConsider)], endChar + 1)

        binLogR <- c()
        for (i in 1:(length(seqToConsider)-1)) {
          PotentialSites <- as.character(seqToConsider[i]:seqToConsider[i+1])
          shared_sites <- PotentialSites %>%
            intersect(names(HLA_A_type1tumorCov)) %>%
            intersect(names(HLA_A_type1normalCov))
          combinedBinTumor <- {
            median(HLA_A_type1tumorCov[shared_sites], na.rm = TRUE) +
            median(HLA_A_type2tumorCov[shared_sites], na.rm = TRUE)
          } %>% repl_NA
          combinedBinNormal <- {
            median(HLA_A_type1normalCov[shared_sites], na.rm= TRUE) +
            median(HLA_A_type2normalCov[shared_sites], na.rm = TRUE)
          } %>% repl_NA
          combinedBinlogR <-
          { log2(combinedBinTumor) - log2(combinedBinNormal) +
              log2(MultFactor) } %>% repl_NA
          type1BinlogR <- shared_sites %>%
            { log2(HLA_A_type1tumorCov[.]) - log2(HLA_A_type1normalCov[.]) +
            log2(MultFactor) } %>%
            median(na.rm = TRUE)
          type2BinlogR <- shared_sites %>%
            { log2(HLA_A_type1tumorCov[.]) - log2(HLA_A_type1normalCov[.]) +
            log2(MultFactor) } %>%
            median(na.rm = TRUE)
          binLogR <- rbind(binLogR,
            cbind(seqToConsider[i], seqToConsider[i+1],
              combinedBinlogR, type1BinlogR, type2BinlogR))
        }

        missMatchseq1_intersect <- as.character(missMatchseq1) %>%
          intersect(names(HLA_A_type1tumorCov)) %>%
          intersect(names(HLA_A_type1normalCov))

        missMatchseq2_intersect <- as.character(missMatchseq2) %>%
          intersect(names(HLA_A_type2tumorCov)) %>%
          intersect(names(HLA_A_type2normalCov))

        ## 2019-04-10 11:10 M.S. was getting warnings in element-wise division
        ## of arrays that were of incompatible sizes. This addition should
        ## prevent that.
        if (length(missMatchseq1_intersect) > 0) {
          type1_logR <-
            element_divide_vector(HLA_A_type1tumorCov[missMatchseq1_intersect],
            HLA_A_type1normalCov[missMatchseq1_intersect],
            mult_factor = MultFactor) %>%
            .[, 2]
        } else {
          type1_logR <- NA
        }

        if (length(missMatchseq2_intersect) > 0) {
          type2_logR <- element_divide_vector(
            HLA_A_type2tumorCov[missMatchseq2_intersect],
            HLA_A_type2normalCov[missMatchseq2_intersect],
            mult_factor = MultFactor) %>%
            set_colnames(c('idx', 'logR_2')) %>%
            .[, 2]
        } else {
          type2_logR <- NA
        }

        tmpOut_cn <- cbind_uneven(
          as.numeric(missMatchseq1_intersect),
          type1_logR,
          HLA_A_type1tumorCov[missMatchseq1_intersect],
          as.numeric(missMatchseq2_intersect),
          type2_logR,
          HLA_A_type2tumorCov[missMatchseq2_intersect],
          HLA_A_type1normalCov[as.character(missMatchseq1_intersect)],
          HLA_A_type2normalCov[as.character(missMatchseq2_intersect)])
        colnames(tmpOut_cn) <-
          c('missMatchseq1', 'logR_type1', 'TumorCov_type1',
            'missMatchseq2', 'logR_type2', 'TumorCov_type2',
            'NormalCov_type1', 'NormalCov_type2')

        dup1 <- unique(tmpOut_cn[duplicated(tmpOut_cn[, 1]), 1])
        dup2 <- unique(tmpOut_cn[duplicated(tmpOut_cn[, 4]), 4])

        ## Average observations in case of duplicates
        for (duplicationIn1 in dup1) {
          tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'TumorCov_type2'] <-
            mean(tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'TumorCov_type2'])
          tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'NormalCov_type2'] <-
            mean(tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'NormalCov_type2'])
          tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'logR_type2'] <-
            mean(tmpOut_cn[tmpOut_cn[, 1] == duplicationIn1, 'logR_type2'])
        }

        for (duplicationIn2 in dup2) {
          tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'TumorCov_type1'] <-
            mean(tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'TumorCov_type1'])
          tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'NormalCov_type1'] <-
            mean(tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'NormalCov_type1'])
          tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'logR_type1'] <-
            mean(tmpOut_cn[tmpOut_cn[, 4] == duplicationIn2, 'logR_type1'])
        }

        tmpOut_cn <- tmpOut_cn[!duplicated(tmpOut_cn[, 1]), , drop = FALSE]
        tmpOut_cn <- tmpOut_cn[!duplicated(tmpOut_cn[, 4]), , drop = FALSE]
        colnames(tmpOut_cn) <-
          c('missMatchseq1', 'logR_type1', 'TumorCov_type1',
            'missMatchseq2', 'logR_type2', 'TumorCov_type2',
            'NormalCov_type1', 'NormalCov_type2')

        combinedTable <- data.frame(tmpOut_cn, stringsAsFactors = FALSE) %>%
          lapply(as.numeric) %>%
          lapply(repl_NA, 0) %>%
          as.data.frame 
        combinedTable$logRcombined <- log2(
          (combinedTable$TumorCov_type1 + combinedTable$TumorCov_type2) /
          (combinedTable$NormalCov_type1 + combinedTable$NormalCov_type2) * 
          MultFactor)
        combinedTable$BAFcombined <- 
          combinedTable$TumorCov_type1 / 
          (combinedTable$TumorCov_type1 + combinedTable$TumorCov_type2)

        if (nrow(combinedTable) != 0) {
          combinedTable$binlogRCombined <- NA
          combinedTable$binlogRtype1 <- NA
          combinedTable$binlogRtype2 <- NA
          combinedTable$binNum <- NA

          ## next, add the binLogR to this table
          for (i in 1:nrow(combinedTable)) {
            ## Retrieve the suitable index in the binLogR table, which contains
            ## LogR estimates obtained using an averaging window
            binned_idx <- {
                (binLogR[, 1] <= combinedTable$missMatchseq1[i]) &
                (binLogR[, 2] > combinedTable$missMatchseq1[i])
              } %>%
              setNames(NULL) %>%
              which

            if (length(binned_idx) == 0 || is.na(binned_idx)) 
              next

            combinedTable[i, ]$binNum <- binLogR[binned_idx, 1]
            combinedTable[i, ]$binlogRCombined <- binLogR[binned_idx, 3]
            combinedTable[i, ]$binlogRtype1 <- binLogR[binned_idx, 4]
            combinedTable[i, ]$binlogRtype2 <- binLogR[binned_idx, 5]
          }
        } else {
          combinedTable$binlogRCombined <- NA
          combinedTable$binlogRtype1 <- NA
          combinedTable$binlogRtype2 <- NA
          combinedTable$binNum <- NA
        }

        rawVals <- funCalcN_withBAF(combinedTable$logRcombined,
          combinedTable$BAFcombined, tumorPloidy, tumorPurity, gamma)
        if (nrow(rawVals) > 0) {
          combinedTable$nAcombined <- rawVals[, 1]
          combinedTable$nBcombined <- rawVals[, 2]
        } else {
          combinedTable$nAcombined <- NA
          combinedTable$nBcombined <- NA
        }

        nB_rawVal_withBAF <- median(combinedTable$nBcombined,
          na.rm = TRUE)
        nB_rawVal_withBAF_conf <- t.test.NA(combinedTable$nBcombined)
        nB_rawVal_withBAF_lower <- nB_rawVal_withBAF_conf$conf.int[1]
        nB_rawVal_withBAF_upper <- nB_rawVal_withBAF_conf$conf.int[2]

        nA_rawVal_withBAF <- median(combinedTable$nAcombined,
          na.rm = TRUE)
        nA_rawVal_withBAF_conf <- t.test.NA(combinedTable$nAcombined)
        nA_rawVal_withBAF_lower <- nA_rawVal_withBAF_conf$conf.int[1]
        nA_rawVal_withBAF_upper <- nA_rawVal_withBAF_conf$conf.int[2]

        rawValsBin <- funCalcN_withBAF(combinedTable$binlogRCombined,
          combinedTable$BAFcombined, tumorPloidy, tumorPurity, gamma)
        combinedTable$nAcombinedBin <- rawValsBin[, 1]
        combinedTable$nBcombinedBin <- rawValsBin[, 2]

        nB_rawVal_withBAF <- median(combinedTable$nBcombined, na.rm = TRUE)
        nB_rawVal_withBAF_conf <- t.test.NA(combinedTable$nBcombined)
        nB_rawVal_withBAF_lower <- nB_rawVal_withBAF_conf$conf.int[1]
        nB_rawVal_withBAF_upper <- nB_rawVal_withBAF_conf$conf.int[2]

        nA_rawVal_withBAF <- median(combinedTable$nAcombined, na.rm = TRUE)
        nA_rawVal_withBAF_conf <- t.test.NA(combinedTable$nAcombined)
        nA_rawVal_withBAF_lower <- nA_rawVal_withBAF_conf$conf.int[1]
        nA_rawVal_withBAF_upper <- nA_rawVal_withBAF_conf$conf.int[2]

        #let's only count non duplicates
        nB_rawVal_withBAF_bin <- median(
          combinedTable[!duplicated(combinedTable$binNum), ]$nBcombinedBin,
          na.rm = TRUE)
        nB_rawVal_withBAF_bin_conf <- t.test.NA(
          combinedTable[!duplicated(combinedTable$binNum), ]$nBcombinedBin)
        nB_rawVal_withBAF_bin_lower <- nB_rawVal_withBAF_bin_conf$conf.int[1]
        nB_rawVal_withBAF_bin_upper <- nB_rawVal_withBAF_bin_conf$conf.int[2]

        nA_rawVal_withBAF_bin <- median(
          combinedTable[!duplicated(combinedTable$binNum), ]$nAcombinedBin,
          na.rm = TRUE)
        nA_rawVal_withBAF_bin_conf <- t.test.NA(
          combinedTable[!duplicated(combinedTable$binNum), ]$nAcombinedBin)
        nA_rawVal_withBAF_bin_lower <- nA_rawVal_withBAF_bin_conf$conf.int[1]
        nA_rawVal_withBAF_bin_upper <- nA_rawVal_withBAF_bin_conf$conf.int[2]

        combinedTable$nAsep <- funCalcN_withoutBAF(
          combinedTable$logR_type1, tumorPloidy, tumorPurity, gamma)
        combinedTable$nAsepBin <- funCalcN_withoutBAF(
          combinedTable$binlogRtype1, tumorPloidy, tumorPurity, gamma)
        combinedTable$nBsep <- funCalcN_withoutBAF(
          combinedTable$logR_type2, tumorPloidy, tumorPurity, gamma)
        combinedTable$nBsepBin <- funCalcN_withoutBAF(
          combinedTable$binlogRtype2, tumorPloidy, tumorPurity, gamma)

        nB_rawVal_withoutBAF <- median(combinedTable$nBsep, na.rm = TRUE)
        nB_rawVal_withoutBAF_conf <- t.test.NA(combinedTable$nBsep)
        nB_rawVal_withoutBAF_lower <- nB_rawVal_withoutBAF_conf$conf.int[1]
        nB_rawVal_withoutBAF_upper <- nB_rawVal_withoutBAF_conf$conf.int[2]

        nA_rawVal_withoutBAF <- median(combinedTable$nAsep, na.rm = TRUE)
        nA_rawVal_withoutBAF_conf <- t.test.NA(combinedTable$nAsep)
        nA_rawVal_withoutBAF_lower <- nA_rawVal_withoutBAF_conf$conf.int[1]
        nA_rawVal_withoutBAF_upper <- nA_rawVal_withoutBAF_conf$conf.int[2]

        nB_rawVal_withoutBAFBin <- median(
          combinedTable[!duplicated(combinedTable$binNum), ]$nBsepBin,
          na.rm = TRUE)
        nB_rawVal_withoutBAFBin_conf <- t.test.NA(
          combinedTable[!duplicated(combinedTable$binNum), ]$nBsepBin)
        nB_rawVal_withoutBAFBin_lower <- nB_rawVal_withoutBAFBin_conf$conf.int[1]
        nB_rawVal_withoutBAFBin_upper <- nB_rawVal_withoutBAFBin_conf$conf.int[2]

        nA_rawVal_withoutBAFBin <- median(
          combinedTable[!duplicated(combinedTable$binNum), ]$nAsepBin,
          na.rm = TRUE)
        nA_rawVal_withoutBAFBin_conf <- t.test.NA(
          combinedTable[!duplicated(combinedTable$binNum), ]$nAsepBin)
        nA_rawVal_withoutBAFBin_lower <-
          nA_rawVal_withoutBAFBin_conf$conf.int[1]
        nA_rawVal_withoutBAFBin_upper <-
          nA_rawVal_withoutBAFBin_conf$conf.int[2]

        ## We can also predict the BAF from our logR
        combinedTable$expectedBAF <-
          (1 - tumorPurity + tumorPurity * combinedTable$nAsep) /
          (2 - 2 * tumorPurity + tumorPurity *
            (combinedTable$nAsep + combinedTable$nBsep))

        # t tests of coverage
        misMatchCoveredInBoth <-
          cbind(ifelse(missMatchPositions$diffSeq1 %in% names(HLA_A_type1normalCov), 1, 0),
            ifelse(missMatchPositions$diffSeq2 %in% names(HLA_A_type2normalCov), 1, 0))
        missMatchseq1 <- missMatchPositions$diffSeq1[rowSums(misMatchCoveredInBoth) == 2]
        missMatchseq2 <- missMatchPositions$diffSeq2[rowSums(misMatchCoveredInBoth) == 2]

        ## We don't want to count the same miss-match multiple times - let's
        ## take an average where this is the case
        type1 <- element_divide_vector(HLA_A_type1tumorCov,
          HLA_A_type1normalCov, missMatchseq1, MultFactor) %>%
          magrittr::set_colnames(c('pos1', 'logR_1'))
        type2 <- element_divide_vector(HLA_A_type2tumorCov,
          HLA_A_type2normalCov, missMatchseq2, MultFactor) %>%
          magrittr::set_colnames(c('pos2', 'logR_2')) %>%
          dplyr::mutate('pos1' = pos2)
        tmpOut <- tryCatch(merge(type1, type2, all = T, by = 'pos1'),
          error = function(e) { print(e); browser() })

        dup1 <- tryCatch(setdiff(unique(tmpOut[duplicated(tmpOut[, 1]), 1]), NA),
          error = function(e) { print(e); browser() })
        dup2 <- tryCatch(setdiff(unique(tmpOut[duplicated(tmpOut[, 3]), 3]), NA),
          error = function(e) { print(e); browser() })

        for (duplicationIn1 in dup1) {
          tryCatch(tmpOut[tmpOut[, 1] == duplicationIn1, 4] <-
            median(tmpOut[tmpOut[, 1] == duplicationIn1, 4]),
          error = function(e) { print(e); browser() })
        }

        for (duplicationIn2 in dup2) {
          tryCatch(
            tmpOut[tmpOut[, 3] == duplicationIn2, 2] <-
              median(tmpOut[tmpOut[, 3] == duplicationIn2, 2])
          , error = function(e) { print(e); browser() })
        }

        tmpOut <- tmpOut[!duplicated(tmpOut[, 1]), , drop = FALSE]
        tmpOut <- tmpOut[!duplicated(tmpOut[, 3]), , drop = FALSE]

        if (nrow(tmpOut) > 1) {
          PairedTtest <- t.test.NA(tmpOut[, 2], tmpOut[, 4], paired = TRUE)
          UnPairedTtest <- t.test.NA(
            remove_NA(tmpOut[, 2]),
            remove_NA(tmpOut[, 4]), paired = FALSE)
        } else {
          PairedTtest <- list(p.value = NA)
          UnPairedTtest <- list(p.value = NA)
        }

        ## t-test of mismatch sites without counting the same read twice
        if (!any(c(length(HLA_A_type1tumorCov_mismatch_unique),
              length(HLA_A_type2tumorCov_mismatch_unique)) == 0)) {
          type1 <- element_divide_vector(HLA_A_type1tumorCov_mismatch_unique,
            HLA_A_type1normalCov_mismatch_unique, missMatchseq1, MultFactor) %>%
            magrittr::set_colnames(c('pos1', 'logR_1'))
          type2 <- element_divide_vector(HLA_A_type2tumorCov_mismatch_unique,
            HLA_A_type2normalCov_mismatch_unique, missMatchseq2, MultFactor) %>%
            magrittr::set_colnames(c('pos2', 'logR_1')) %>%
            dplyr::mutate('pos1' = pos2)
          tmpOut_unique <- tryCatch(merge(type1, type2, all = T, by = 'pos1'),
          error = function(e) { print(e); browser() })

          dup1_unique <- setdiff(unique(tmpOut_unique[duplicated(tmpOut_unique[,
                1]), 1]), NA)
          dup2_unique <- setdiff(unique(tmpOut_unique[duplicated(tmpOut_unique[,
                3]), 3]), NA)

          for (duplicationIn1 in dup1_unique) {
            tmpOut_unique[tmpOut_unique[, 1] == duplicationIn1, 4] <-
              median(tmpOut_unique[tmpOut_unique[, 1] == duplicationIn1, 4])
          }

          for (duplicationIn2 in dup2_unique) {
            tmpOut_unique[tmpOut_unique[, 3] == duplicationIn2, 2] <-
              median(tmpOut_unique[tmpOut_unique[, 3] == duplicationIn2, 2])
          }

          tmpOut_unique <- tmpOut_unique[!duplicated(tmpOut_unique[, 1]), ,
            drop = FALSE]
          tmpOut_unique <- tmpOut_unique[!duplicated(tmpOut_unique[, 3]), ,
            drop = FALSE]

          if (nrow(tmpOut_unique) > 1) {
            PairedTtest_unique <- t.test.NA(tmpOut_unique[, 2],
              tmpOut_unique[, 4], paired = TRUE)
            UnPairedTtest_unique <- t.test.NA(
              remove_NA(tmpOut_unique[, 2]),
              remove_NA(tmpOut_unique[, 4]), paired = FALSE)
          } else {
            PairedTtest_unique <- list(p.value = NA)
            UnPairedTtest_unique <- list(p.value = NA)
          }
        } else {
          PairedTtest_unique <- list(p.value = NA)
          UnPairedTtest_unique <- list(p.value = NA)
        }

        ## Let's put togehter the output
        HLAtype1Log2MedianCoverage <-
          element_divide_vector(HLA_A_type1tumorCov, HLA_A_type1normalCov,
            mult_factor = MultFactor)[, 2] %>% median(na.rm = T)
        HLAtype2Log2MedianCoverage <-
          element_divide_vector(HLA_A_type2tumorCov, HLA_A_type2normalCov,
            mult_factor = MultFactor)[, 2] %>% median(na.rm = T)
        HLAtype1Log2MedianCoverageAtSites <-
          element_divide_vector(HLA_A_type1tumorCov, HLA_A_type1normalCov,
            sites = missMatchseq1,
            mult_factor = MultFactor)[, 2] %>% median(na.rm = T)
        HLAtype2Log2MedianCoverageAtSites <-
          element_divide_vector(HLA_A_type2tumorCov, HLA_A_type2normalCov,
            sites = missMatchseq2,
            mult_factor = MultFactor)[, 2] %>% median(na.rm = T)
        PVal <- PairedTtest$p.value
        UnPairedPval <- UnPairedTtest$p.value
        PVal_unique <- PairedTtest_unique$p.value
        UnPairedPval_unique <- UnPairedTtest_unique$p.value
        LossAllele <- ifelse(HLAtype1Log2MedianCoverageAtSites <
            HLAtype2Log2MedianCoverageAtSites, 1, 2) %>%
          { c(HLA_A_type1, HLA_A_type2)[.] }
        KeptAllele <- setdiff(c(HLA_A_type1, HLA_A_type2), LossAllele)
        numMisMatchSitesCov <- which(rowSums(misMatchCoveredInBoth) == 2) %>%
          length
        propSupportiveSites <-
          tryCatch(
            {
              a <- element_divide_vector(HLA_A_type1tumorCov, HLA_A_type1normalCov,
                       missMatchseq1)
              b <- element_divide_vector(HLA_A_type2tumorCov, HLA_A_type2normalCov,
                       missMatchseq2)
              min_length <- min(nrow(a), nrow(b), na.rm = T)
              unlist(a[1:min_length, 2]) - unlist(b[1:min_length, 2])
            },
            error = function(e) { print(e); browser() },
            warning = function(e) { print(e); browser() }) %>%
          setdiff(NaN) %>% { . > 0 } %>% mean %>% { max(., 1 - .) * 100 }

        # Additional features we can add in relation to integer copy numbers
        if (performIntegerCopyNum) {
          HLA_type1copyNum_withoutBAF <- nA_rawVal_withoutBAF
          HLA_type1copyNum_withoutBAF_lower <- nA_rawVal_withoutBAF_lower
          HLA_type1copyNum_withoutBAF_upper <- nA_rawVal_withoutBAF_upper

          HLA_type1copyNum_withBAF <- nA_rawVal_withBAF
          HLA_type1copyNum_withBAF_lower <- nA_rawVal_withBAF_lower
          HLA_type1copyNum_withBAF_upper <- nA_rawVal_withBAF_upper

          HLA_type2copyNum_withoutBAF <- nB_rawVal_withoutBAF
          HLA_type2copyNum_withoutBAF_lower <- nB_rawVal_withoutBAF_lower
          HLA_type2copyNum_withoutBAF_upper <- nB_rawVal_withoutBAF_upper

          HLA_type2copyNum_withBAF <- nB_rawVal_withBAF
          HLA_type2copyNum_withBAF_lower <- nB_rawVal_withBAF_lower
          HLA_type2copyNum_withBAF_upper <- nB_rawVal_withBAF_upper

          HLA_type1copyNum_withoutBAFBin <- nA_rawVal_withoutBAFBin
          HLA_type1copyNum_withoutBAFBin_lower <- nA_rawVal_withoutBAFBin_lower
          HLA_type1copyNum_withoutBAFBin_upper <- nA_rawVal_withoutBAFBin_upper

          HLA_type1copyNum_withBAFBin <- nA_rawVal_withBAF_bin
          HLA_type1copyNum_withBAFBin_lower <- nA_rawVal_withBAF_bin_lower
          HLA_type1copyNum_withBAFBin_upper <- nA_rawVal_withBAF_bin_upper

          HLA_type2copyNum_withoutBAFBin <- nB_rawVal_withoutBAFBin
          HLA_type2copyNum_withoutBAFBin_lower <- nB_rawVal_withoutBAFBin_lower
          HLA_type2copyNum_withoutBAFBin_upper <- nB_rawVal_withoutBAFBin_upper

          HLA_type2copyNum_withBAFBin <- nB_rawVal_withBAF_bin
          HLA_type2copyNum_withBAFBin_lower <- nB_rawVal_withBAF_bin_lower
          HLA_type2copyNum_withBAFBin_upper <- nB_rawVal_withBAF_bin_upper
        } else {
          HLA_type1copyNum_withoutBAF <- NA
          HLA_type1copyNum_withoutBAF_lower <- NA
          HLA_type1copyNum_withoutBAF_upper <- NA

          HLA_type1copyNum_withBAF <- NA
          HLA_type1copyNum_withBAF_lower <- NA
          HLA_type1copyNum_withBAF_upper <- NA

          HLA_type2copyNum_withoutBAF <- NA
          HLA_type2copyNum_withoutBAF_lower <- NA
          HLA_type2copyNum_withoutBAF_upper <- NA

          HLA_type2copyNum_withBAF <- NA
          HLA_type2copyNum_withBAF_lower <- NA
          HLA_type2copyNum_withBAF_upper <- NA

          HLA_type1copyNum_withoutBAFBin <- NA
          HLA_type1copyNum_withoutBAFBin_lower <- NA
          HLA_type1copyNum_withoutBAFBin_upper <- NA

          HLA_type1copyNum_withBAFBin <- NA
          HLA_type1copyNum_withBAFBin_lower <- NA
          HLA_type1copyNum_withBAFBin_upper <- NA

          HLA_type2copyNum_withoutBAFBin <- NA
          HLA_type2copyNum_withoutBAFBin_lower <- NA
          HLA_type2copyNum_withoutBAFBin_upper <- NA

          HLA_type2copyNum_withBAFBin <- NA
          HLA_type2copyNum_withBAFBin_lower <- NA
          HLA_type2copyNum_withBAFBin_upper <- NA
        }


        ## save some temporary files before plotting
        statistics_fn <- paste(figureDir, '/', sample, '.', HLA_gene,
            '.tmp.data.plots.RData', sep = '')
        save.image(statistics_fn)

        return(list(
          message='analysis_ok',
          HLA_A_type1=HLA_A_type1,
          HLA_A_type2=HLA_A_type2,
          HLAtype1Log2MedianCoverage=HLAtype1Log2MedianCoverage,
          HLAtype2Log2MedianCoverage=HLAtype2Log2MedianCoverage,
          HLAtype1Log2MedianCoverageAtSites=HLAtype1Log2MedianCoverageAtSites,
          HLAtype2Log2MedianCoverageAtSites=HLAtype2Log2MedianCoverageAtSites,
          HLA_type1copyNum_withoutBAF=HLA_type1copyNum_withoutBAF,
          HLA_type1copyNum_withoutBAF_lower=HLA_type1copyNum_withoutBAF_lower,
          HLA_type1copyNum_withoutBAF_upper=HLA_type1copyNum_withoutBAF_upper,
          HLA_type1copyNum_withBAF=HLA_type1copyNum_withBAF,
          HLA_type1copyNum_withBAF_lower=HLA_type1copyNum_withBAF_lower,
          HLA_type1copyNum_withBAF_upper=HLA_type1copyNum_withBAF_upper,
          HLA_type2copyNum_withoutBAF=HLA_type2copyNum_withoutBAF,
          HLA_type2copyNum_withoutBAF_lower=HLA_type2copyNum_withoutBAF_lower,
          HLA_type2copyNum_withoutBAF_upper=HLA_type2copyNum_withoutBAF_upper,
          HLA_type2copyNum_withBAF=HLA_type2copyNum_withBAF,
          HLA_type2copyNum_withBAF_lower=HLA_type2copyNum_withBAF_lower,
          HLA_type2copyNum_withBAF_upper=HLA_type2copyNum_withBAF_upper,
          HLA_type1copyNum_withoutBAFBin=HLA_type1copyNum_withoutBAFBin,
          HLA_type1copyNum_withoutBAFBin_lower=HLA_type1copyNum_withoutBAFBin_lower,
          HLA_type1copyNum_withoutBAFBin_upper=HLA_type1copyNum_withoutBAFBin_upper,
          HLA_type1copyNum_withBAFBin=HLA_type1copyNum_withBAFBin,
          HLA_type1copyNum_withBAFBin_lower=HLA_type1copyNum_withBAFBin_lower,
          HLA_type1copyNum_withBAFBin_upper=HLA_type1copyNum_withBAFBin_upper,
          HLA_type2copyNum_withoutBAFBin=HLA_type2copyNum_withoutBAFBin,
          HLA_type2copyNum_withoutBAFBin_lower=HLA_type2copyNum_withoutBAFBin_lower,
          HLA_type2copyNum_withoutBAFBin_upper=HLA_type2copyNum_withoutBAFBin_upper,
          HLA_type2copyNum_withBAFBin=HLA_type2copyNum_withBAFBin,
          HLA_type2copyNum_withBAFBin_lower=HLA_type2copyNum_withBAFBin_lower,
          HLA_type2copyNum_withBAFBin_upper=HLA_type2copyNum_withBAFBin_upper,
          PVal=PVal,
          UnPairedPval=UnPairedPval,
          PVal_unique=PVal_unique,
          UnPairedPval_unique=UnPairedPval_unique,
          LossAllele=LossAllele,
          KeptAllele=KeptAllele,
          numMisMatchSitesCov=numMisMatchSitesCov,
          combinedTable=combinedTable,
          propSupportiveSites=propSupportiveSites))
      }
    })

    HLAoutPut <- tryCatch(
      plyr::llply(HLAoutPut_l, function(x) {
        x <- x[setdiff(names(x), 'combinedTable')]
        x[sapply(x, is.null)] <- NA;
        class(x$message) <- 'character'
        x
      }) %>% rbindlist(fill = T)
      , error = function(e) { 
      print('Problem with HLAoutPut'); print(e); browser() 
    })

    combinedTable <- tryCatch(purrr::imap(HLAoutPut_l, function(x, idx)
        as.data.frame(cbind('hla' = hlas[idx], x[['combinedTable']]))), 
      error = function(e) { 
        print('Problem with combinedTable'); print(e); NULL 
      })

    combinedTable <- tryCatch(rbindlist(combinedTable, fill = T),
      error = function(e) { 
        print('Problem with combinedTable'); print(e); NULL 
      })
  }
  ### Coverage step }}}

  ## {{{ Plotting step
  if (plottingStep) {
    logger('plotting for sample: ', sample)
    min_coverage_fn <-
      paste(figureDir, '/', sample, '.minCoverage_', minCoverageFilter,
        '.HLA.pdf', sep = '')
    pdf(min_coverage_fn, width = 10, height = 6)

    for (HLA_gene in unique(substr(hlaAlleles, 1, 5))) {
      # if (file.exists(paste(figureDir, '/', sample, '.',
      #       HLA_gene, '.tmp.data.plots.RData', sep = ''))) {
      #   load(paste(figureDir, '/', sample, '.', HLA_gene,
      #       '.tmp.data.plots.RData', sep = ''))
      # }

      if (F) {
        msg <- paste('Run with coverageStep == TRUE for : ', sample,
          ' first!', sep = '')
        howToWarn(msg)
        logger(msg)
        next
      }

      # if (!exists('regionSpecOutPut')) {
      #   msg <- paste('\ncoverageStep did not run to completion for: ',
      #     HLA_gene, ' in ', sample, '! ', '\n', sep = '')
      #   logger(msg)
      #   stop(msg)
      # }

      ## getting exons for both alleles
      ## {{{ Type 1
      HLA_A_type1DatFormat <- unlist(strsplit(HLA_A_type1, split = "_")) %>%
        .[2:length(unlist(strsplit(HLA_A_type1, split = "_")))]
      HLA_A_type1Formatted <- toupper(paste('HLA-', HLA_A_type1DatFormat[1],
          '\\*',
          paste(HLA_A_type1DatFormat[2:length(HLA_A_type1DatFormat)],
            collapse = ':'), sep = ''))
      cmd <- paste('awk \'/^DE   ', HLA_A_type1Formatted,
        ", / {p = 1}; p; /Sequence/ {p = 0}\' ", HLAexonLoc, sep = '')
      awk.result <- system(cmd, intern = TRUE)
      HLAtype1exons <- grep('^FT   exon', awk.result, value = TRUE)
      HLAtype1exons_s <- strsplit(HLAtype1exons, split = " ")
      HLAtype1exons <- do.call(rbind, HLAtype1exons_s)
      HLAtype1exons <- HLAtype1exons[, ncol(HLAtype1exons)]
      if (length(HLAtype1exons) != 0) {
        HLAtype1exonTable <- c()
        for (i in 1:length(HLAtype1exons)) {
          HLAtype1exonTable <- rbind(HLAtype1exonTable,
            (unlist(strsplit(HLAtype1exons[i], split = "\\.\\."))))
        }
      }
      ## }}} Type 1
      ## {{{ Type 2
      HLA_A_type2DatFormat <- unlist(strsplit(HLA_A_type2, split = "_")) %>%
        .[2:length(unlist(strsplit(HLA_A_type2, split = "_")))]
      HLA_A_type2Formatted <- toupper(paste('HLA-', HLA_A_type2DatFormat[2],
          '\\*',
          paste(HLA_A_type2DatFormat[2:length(HLA_A_type2DatFormat)],
            collapse = ':'), sep = ''))
      cmd <- paste("awk \'/^DE   ", HLA_A_type2Formatted,
        ", / {p = 2}; p; /Sequence/ {p = 0}\' ", HLAexonLoc, sep = '')
      awk.result <- system(cmd, intern = TRUE)
      HLAtype2exons <- grep('^FT   exon', awk.result, value = TRUE)
      HLAtype2exons_s <- strsplit(HLAtype2exons, split = " ")
      HLAtype2exons <- do.call(rbind, HLAtype2exons_s)
      HLAtype2exons <- HLAtype2exons[, ncol(HLAtype2exons)]
      if (length(HLAtype2exons) != 0) {
        HLAtype2exonTable <- c()
        for (i in 2:length(HLAtype2exons)) {
          HLAtype2exonTable <- rbind(HLAtype2exonTable,
            (unlist(strsplit(HLAtype2exons[i], split = "\\.\\."))))
        }
      }
      ## }}} Type 2

      ## Some things to plot if there is a normal sample
      if (runWithNormal) {
        ## Rolling mean
        tryCatch({
        par(mfrow = c(2, 1))
        par(mar = c(2, 5, 2, 2))
        barplot(c(rollmean(HLA_A_type2tumorCov*MultFactor, min(500, length(HLA_A_type2tumorCov)))/rollmean(HLA_A_type2normalCov, min(500, length(HLA_A_type2normalCov)))), ylim = c(0, 3), xaxt = 'n', main = HLA_A_type2, las = 1, ylab = 'Tumor/Normal Coverage')
        abline(h = 1, lty = 'dashed', col = 'blue', lwd = 1.5)
        barplot(c(rollmean(HLA_A_type1tumorCov*MultFactor,
          min(500, length(HLA_A_type1tumorCov))) /
            rollmean(HLA_A_type1normalCov,
          min(500, length(HLA_A_type1normalCov)))),
          ylim = c(0, 3), xaxt = 'n', main = HLA_A_type1, las = 1,
          ylab = 'Tumor/Normal Coverage')
        abline(h = 1, lty = 'dashed', col = 'blue', lwd = 1.5)
        }, error = function(e) return(NULL))

        ## Log ratio and density of mismatches
        tryCatch({
        par(mfrow = c(1, 1))
        par(mar = c(5, 5, 5, 2))
        plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2,
              HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)),
          lim = c(-3, 3), col = '#3182bd99', pch = 16,
          lab = 'HLA genomic position',
          lab = 'Log Ratio',
          ain = c(paste("HLA raw balance", sample)),
          ex = 0.75, ype = 'n' , as = 1)
        ## add the exonic positions
        if (length(HLAtype1exons) != 0) {
          for (i in 1:nrow(HLAtype1exonTable)) {
            rect(xleft = as.numeric(HLAtype1exonTable[i, 1]),
              xright = as.numeric(HLAtype1exonTable[i, 2]),
              ybottom = -2, ytop = 2, col = '#bdbdbd25', border = FALSE)
          }
        }

        if (length(HLAtype2exons) != 0) {
          for (i in 1:nrow(HLAtype2exonTable)) {
            rect(xleft = as.numeric(HLAtype2exonTable[i, 1]),
              xright = as.numeric(HLAtype2exonTable[i, 2]),
              ybottom = -2, ytop = 2, col = '#bdbdbd25', border = FALSE)
          }
        }

        points(c(names(HLA_A_type2normalCov)),
          log2(c((HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)),
          col = '#3182bd99', pch = 16)
        points(names(HLA_A_type2normalCov)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], col = 'black', bg = '#3182bd99', pch = 21, cex = 1)
        points(names(HLA_A_type1normalCov), log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor), col = '#de2d2699', pch = 16, cex = 0.75)
        points(names(HLA_A_type1normalCov)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], col = 'black', bg = '#de2d2699', pch = 21, cex = 1)

        points(missMatchPositions$diffSeq1,
          rep(-3, length(missMatchPositions$diffSeq1)),
          col = '#63636399', bg = '#63636399', pch = 21, cex = 1)
        d <- density(missMatchPositions$diffSeq1, bw = 40)
        d$y <- (d$y/max(d$y))
        d$y <- d$y -3
        lines(d)
        abline(h = 0, lty = 'dashed')

        legend('topright', legend = c(HLA_A_type2, HLA_A_type1),
               lty = 1, col = c('#3182bd99', '#de2d2699'), bty = 'n',
               cex = 1, lwd = 3)
        }, error = function(e) return(NULL))

        tryCatch({
        # normal coverage comparison
        plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2))
             , lim = c(0, max(HLA_A_type1normalCov, HLA_A_type2normalCov)), col = '#3182bd99', pch = 16
             , lab = 'HLA genomic position'
             , lab = 'Coverage'
             , ain = c(paste("HLA normal coverage", sample))
             , ex = 0.75
             , ype = 'n'
             , as = 1)
        # add the exonic positions
        if (length(HLAtype1exons)!= 0) {
          for (i in 1:nrow(HLAtype1exonTable)) {
            rect(xleft = as.numeric(HLAtype1exonTable[i, 1]), xright = as.numeric(HLAtype1exonTable[i, 2]), ybottom = -2, ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)), col = '#bdbdbd25', border = FALSE)
          }
        }

        if (length(HLAtype2exons)!= 0) {
          for (i in 1:nrow(HLAtype2exonTable)) {
            rect(xleft = as.numeric(HLAtype2exonTable[i, 1]), xright = as.numeric(HLAtype2exonTable[i, 2]), ybottom = -2, ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)), col = '#bdbdbd25', border = FALSE)
          }
        }

        lines(c(HLA_A_type2normal$V2), c(HLA_A_type2normalCov), col = '#3182bd')
        lines(c(HLA_A_type1normal$V2), c(HLA_A_type1normalCov), col = '#de2d26')
        legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
               lty = 1, col = c('#3182bd', '#de2d26'), bty = 'n', cex = 1, lwd = 3)

        points(c(HLA_A_type1normal$V2)[HLA_A_type1normal$V2 %in% missMatchPositions$diffSeq1], c(HLA_A_type1normalCov)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], col = '#de2d26', pch = 16)
        points(c(HLA_A_type2normal$V2)[HLA_A_type2normal$V2 %in% missMatchPositions$diffSeq2], c(HLA_A_type2normalCov)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], col = '#3182bd', pch = 16)
        }, error = function(e) return(NULL))

        tryCatch({
        # tumor and normal coverage for allele 1
        plot(names(HLA_A_type1normalCov), apply(cbind(HLA_A_type1tumorCov*MultFactor, HLA_A_type1normalCov), 1, max), type = 'n', xlab = 'HLA genomic position', ylab = 'Coverage', las = 1, main = HLA_A_type1, ylim = c(0, max(c(HLA_A_type1tumorCov, HLA_A_type1normalCov))))
        if (length(HLAtype1exons)!= 0)
        {
          for (i in 1:nrow(HLAtype1exonTable))
          {
            rect(xleft = as.numeric(HLAtype1exonTable[i, 1]), xright = as.numeric(HLAtype1exonTable[i, 2]), ybottom = 0, ytop = max(c(HLA_A_type1tumorCov, HLA_A_type1normalCov)), col = '#bdbdbd50', border = FALSE)
          }
        }

        lines(c(HLA_A_type1normal$V2), c(HLA_A_type1tumorCov*MultFactor), col = '#3182bd')
        lines(c(HLA_A_type1normal$V2), c(HLA_A_type1normalCov), col = '#9ecae1')
        legend('topleft', legend = c('tumour', 'normal') ,
               lty = 1, col = c('#3182bd', '#9ecae1'), bty = 'n', cex = 1, lwd = 3)

        points(c(HLA_A_type1normal$V2)[HLA_A_type1normal$V2 %in% missMatchPositions$diffSeq1], c(HLA_A_type1tumorCov*MultFactor)[names(HLA_A_type1tumorCov) %in% missMatchPositions$diffSeq1], col = '#3182bd', pch = 16)
        }, error = function(e) return(NULL))

        tryCatch({
        # tumor and normal coverage for allele 2
        plot(c(HLA_A_type2normal$V2), apply(cbind(HLA_A_type2tumorCov*MultFactor, HLA_A_type2normalCov), 1, max), type = 'n', xlab = 'HLA genomic position', ylab = 'Coverage', las = 1, main = HLA_A_type2, ylim = c(0, max(c(HLA_A_type2tumorCov, HLA_A_type2normalCov))))
        if (length(HLAtype2exons)!= 0)
        {
          for (i in 1:nrow(HLAtype2exonTable))
          {
            rect(xleft = as.numeric(HLAtype2exonTable[i, 1]), xright = as.numeric(HLAtype2exonTable[i, 2]), ybottom = 0, ytop = max(c(HLA_A_type2tumorCov, HLA_A_type2normalCov)), col = '#bdbdbd50', border = FALSE)
          }
        }

        lines(c(HLA_A_type2normal$V2), c(HLA_A_type2tumorCov*MultFactor), col = '#3182bd')
        lines(c(HLA_A_type2normal$V2), c(HLA_A_type2normalCov), col = '#9ecae1')
        legend('topleft', legend = c('tumour', 'normal') ,
               lty = 1, col = c('#3182bd', '#9ecae1'), bty = 'n', cex = 1, lwd = 3)
        # points(c(HLA_A_type2normal$V2)[HLA_A_type2tumor$V2 %in% missMatchPositions$diffSeq2], c(HLA_A_type2tumorCov)[names(HLA_A_type2tumorCov) %in% missMatchPositions$diffSeq2], col = '#3182bd', pch = 16)
        points(c(HLA_A_type2normal$V2)[HLA_A_type2normal$V2 %in% missMatchPositions$diffSeq2], c(HLA_A_type2tumorCov*MultFactor)[names(HLA_A_type2tumorCov*MultFactor) %in% missMatchPositions$diffSeq2], col = '#3182bd', pch = 16)
        }, error = function(e) return(NULL))



        tryCatch({
        #let's now just plot the mismatch positions
        plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2, HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)), ylim = c(-2, 2), col = '#3182bd99', pch = 16
             , lab = 'HLA genomic position'
             , lab = 'Log Ratio'
             , ain = c(paste("HLA raw balance", sample))
             , ex = 0.75
             , ype = 'n')
        points(c(HLA_A_type2normal$V2)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2],
          log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], col = '#3182bd99', pch = 16, cex = 1)
        points(c(HLA_A_type1normal$V2)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1],
          log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], col = '#de2d2699', pch = 16, cex = 1)
        abline(h = 0, lty = 'dashed')
        }, error = function(e) return(NULL))

      }

      ## fewer things to plot if there's not a normal sample
      if (!runWithNormal) {
        tryCatch({
        par(mfrow = c(2, 1))
        par(mar = c(2, 5, 2, 2))
        Ymax <- max(c(HLA_A_type1tumorCov, HLA_A_type2tumorCov))+100

        barplot(c(rollmean(HLA_A_type1tumorCov, 150)),
          ylim = c(0, Ymax), xaxt = 'n', main = HLA_A_type1,
          las = 1, col = '#de2d2699', border = '#de2d2650')
        abline(h = median(HLA_A_type1tumorCov), lty = 'dashed',
          col = '#de2d26', lwd = 1.5, na.rm = TRUE)
        abline(h = median(HLA_A_type2tumorCov), lty = 'dashed',
          col = '#3182bd', lwd = 1.5, na.rm = TRUE)
        barplot(c(rollmean(HLA_A_type2tumorCov, 150)),
          ylim = c(0, Ymax), xaxt = 'n', main = HLA_A_type2,
          las = 1, col = '#3182bd', border = '#3182bd50')
        abline(h = median(HLA_A_type1tumorCov), lty = 'dashed',
          col = '#de2d26', lwd = 1.5, na.rm = TRUE)
        abline(h = median(HLA_A_type2tumorCov), lty = 'dashed',
          col = '#3182bd', lwd = 1.5, na.rm = TRUE)
        }, error = function(e) return(NULL))

        tryCatch({
        par(mfrow = c(1, 1))
        par(mar = c(5, 5, 5, 2))
        plot(c(1:max(c(HLA_A_type1tumor$V2, HLA_A_type2tumor$V2)))
             , lim = c(0, Ymax), col = '#3182bd99', pch = 16
             , lab = 'HLA genomic position'
             , lab = 'Coverage'
             , ain = c(paste("HLA raw coverage", sample))
             , ex = 0.75
             , ype = 'n'
             , as = 1)

        points(HLA_A_type1tumor$V2, HLA_A_type1tumor$V4, col = '#de2d2699', pch = 16)
        points(HLA_A_type1tumor$V2[HLA_A_type1tumor$V2 %in% missMatchPositions$diffSeq1], HLA_A_type1tumor$V4[HLA_A_type1tumor$V2 %in% missMatchPositions$diffSeq1], col = 'black', bg = '#de2d2699', pch = 21, cex = 1)

        points(HLA_A_type2tumor$V2, HLA_A_type2tumor$V4, col = '#3182bd99', pch = 16)
        points(HLA_A_type2tumor$V2[HLA_A_type2tumor$V2 %in% missMatchPositions$diffSeq2], HLA_A_type2tumor$V4[HLA_A_type2tumor$V2 %in% missMatchPositions$diffSeq2], col = 'black', bg = '#3182bd99', pch = 21, cex = 1)

        legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
               lty = 1, col = c('#3182bd99', '#de2d2699'), bty = 'n', cex = 1, lwd = 3)
        }, error = function(e) return(NULL))


        tryCatch({
        plot(c(1:max(c(HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)))
             , lim = c(0, Ymax), col = '#3182bd99', pch = 16
             , lab = 'HLA genomic position'
             , lab = 'Coverage'
             , ain = c(paste("HLA raw balance", sample))
             , ex = 0.75
             , ype = 'n'
             , as = 1)
        points(HLA_A_type1tumor$V2[HLA_A_type1tumor$V2 %in% missMatchPositions$diffSeq1], HLA_A_type1tumor$V4[HLA_A_type1tumor$V2 %in% missMatchPositions$diffSeq1], col = '#de2d2699', pch = 16, cex = 1)
        points(HLA_A_type2tumor$V2[HLA_A_type2tumor$V2 %in% missMatchPositions$diffSeq2], HLA_A_type2tumor$V4[HLA_A_type2tumor$V2 %in% missMatchPositions$diffSeq2], col = '#3182bd99', pch = 16, cex = 1)

        legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
               lty = 1, col = c('#3182bd99', '#de2d2699'), bty = 'n', cex = 1, lwd = 3)
        }, error = function(e) return(NULL))

      }

      ## next, let's plot the predicted integer copy numbers
      if (!is.na(tumorPurity) && nrow(combinedTable) != 0) {
        tryCatch({
        plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2, HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)), ylim = c(-0.5, max(round(c(combinedTable$nBcombined, combinedTable$nAcombined)), na.rm = TRUE)+1), col = '#3182bd99', pch = 16
             , lab = 'HLA genomic position'
             , lab = 'Copy Number'
             , ain = c(paste("HLA copyNum balance", sample))
             , ex = 0.75
             , ype = 'n'
             , axt = 'n')

        if (length(HLAtype1exons)!= 0) {
          for (i in 1:nrow(HLAtype1exonTable)) {
            rect(xleft = as.numeric(HLAtype1exonTable[i, 1]), xright = as.numeric(HLAtype1exonTable[i, 2]), ybottom = -2, ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)), col = '#bdbdbd25', border = FALSE)
          }
        }

        if (length(HLAtype2exons)!= 0) {
          for (i in 1:nrow(HLAtype2exonTable)) {
            rect(xleft = as.numeric(HLAtype2exonTable[i, 1]), xright = as.numeric(HLAtype2exonTable[i, 2]), ybottom = -2, ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)), col = '#bdbdbd25', border = FALSE)
          }
        }

        if (useLogRbin) {
          axis(side = 2, at = 0:max(round(c(combinedTable$nAcombinedBin, combinedTable$nBcombinedBin)), na.rm = TRUE), las = 1)
          points(combinedTable$missMatchseq1, combinedTable$nAcombinedBin, col = '#de2d2699', pch = 16, cex = 1)
          points(combinedTable$missMatchseq2, combinedTable$nBcombinedBin, col = '#3182bd99', pch = 16, cex = 1)

          abline(h = nA_rawVal_withBAF_bin, lty = 'dashed', col = '#de2d2699', lwd = 3)
          abline(h = nB_rawVal_withBAF_bin, lty = 'dashed', col = '#3182bd99', lwd = 3)

          legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
                 lty = 1, col = c('#3182bd', '#de2d26'), bty = 'n', cex = 1, lwd = 3)
        }

        if (!useLogRbin) {
          axis(side = 2, at = 0:max(round(c(combinedTable$nAcombined, combinedTable$nBcombined)), na.rm = TRUE), las = 1)
          points(combinedTable$missMatchseq1, combinedTable$nAcombined, col = '#de2d2699', pch = 16, cex = 1)
          points(combinedTable$missMatchseq2, combinedTable$nBcombined, col = '#3182bd99', pch = 16, cex = 1)

          abline(h = nA_rawVal_withBAF, lty = 'dashed', col = '#de2d2699', lwd = 3)
          abline(h = nB_rawVal_withBAF, lty = 'dashed', col = '#3182bd99', lwd = 3)

          legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
                 lty = 1, col = c('#3182bd', '#de2d26'), bty = 'n', cex = 1, lwd = 3)
        }
        }, error = function(e) return(NULL))

      }

      if (!is.na(tumorPurity) && nrow(combinedTable) != 0) {
        if (!useLogRbin) {
          tryCatch({
          plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2,
                HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)),
            ylim = c(-0.5, max(round(c(combinedTable$nAsep,
                    combinedTable$nBsep)), na.rm = TRUE) + 1),
            col = '#3182bd99', pch = 16, lab = 'HLA genomic position',
            lab = 'Copy Number',
            ain = c(paste("HLA copyNum balance", sample)),
            ex = 0.75, ype = 'n', axt = 'n')
          axis(side = 2,
            at = 0:max(round(c(combinedTable$nAsep, combinedTable$nBsep))),
            las = 1)

          if (length(HLAtype1exons)!= 0) {
            for (i in 1:nrow(HLAtype1exonTable)) {
              rect(xleft = as.numeric(HLAtype1exonTable[i, 1]),
                xright = as.numeric(HLAtype1exonTable[i, 2]),
                ybottom = -2,
                ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),
                col = '#bdbdbd25', border = FALSE)
            }
          }

          if (length(HLAtype2exons)!= 0) {
            for (i in 1:nrow(HLAtype2exonTable)) {
              rect(xleft = as.numeric(HLAtype2exonTable[i, 1]),
                xright = as.numeric(HLAtype2exonTable[i, 2]),
                ybottom = -2,
                ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),
                col = '#bdbdbd25', border = FALSE)
            }
          }

          points(combinedTable$missMatchseq1, combinedTable$nAsep,
            col = '#de2d2699', pch = 16, cex = 1)
          points(combinedTable$missMatchseq2, combinedTable$nBsep,
            col = '#3182bd99', pch = 16, cex = 1)

          abline(h = nA_rawVal_withoutBAF, lty = 'dashed',
            col = '#de2d2699', lwd = 3)
          abline(h = nB_rawVal_withoutBAF, lty = 'dashed',
            col = '#3182bd99', lwd = 3)

          legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
            lty = 1, col = c('#3182bd', '#de2d26'),
            bty = 'n', cex = 1, lwd = 3)
          }, error = function(e) return(NULL))
        }

        if (useLogRbin) {
          tryCatch({
          plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2,
                HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)),
            ylim = c(-0.5, max(round(c(combinedTable$nAsep, combinedTable$nBsep)),
                na.rm = TRUE) + 1), col = '#3182bd99', pch = 16,
            lab = 'HLA genomic position', lab = 'Copy Number',
            ain = c(paste("HLA copyNum balance", sample)),
            ex = 0.75, ype = 'n', axt = 'n')
          axis(side = 2,
            at = 0:max(round(c(combinedTable$nAsep, combinedTable$nBsep))),
            las = 1)

          if (length(HLAtype1exons) != 0) {
            for (i in 1:nrow(HLAtype1exonTable)) {
              rect(xleft = as.numeric(HLAtype1exonTable[i, 1]),
                xright = as.numeric(HLAtype1exonTable[i, 2]),
                ybottom = -2,
                ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),
                col = '#bdbdbd25', border = FALSE)
            }
          }

          if (length(HLAtype2exons) != 0) {
            for (i in 1:nrow(HLAtype2exonTable)) {
              rect(xleft = as.numeric(HLAtype2exonTable[i, 1]),
                xright = as.numeric(HLAtype2exonTable[i, 2]),
                ybottom = -2,
                ytop = max(c(HLA_A_type1normalCov, HLA_A_type2normalCov)),
                col = '#bdbdbd25', border = FALSE)
            }
          }

          points(combinedTable$missMatchseq1,
            combinedTable$nAsepBin, col = '#de2d2699', pch = 16, cex = 1)
          points(combinedTable$missMatchseq2,
            combinedTable$nBsepBin, col = '#3182bd99', pch = 16, cex = 1)

          abline(h = nA_rawVal_withoutBAFBin, lty = 'dashed',
            col = '#de2d2699', lwd = 3)
          abline(h = nB_rawVal_withoutBAFBin, lty = 'dashed',
            col = '#3182bd99', lwd = 3)

          legend('topleft', legend = c(HLA_A_type2, HLA_A_type1) ,
            lty = 1, col = c('#3182bd', '#de2d26'),
            bty = 'n', cex = 1, lwd = 3)
          }, error = function(e) return(NULL))
        }
      }

      # t-test plot
      par(mfrow = c(1, 2))
      par(mar = c(5, 5, 5, 2))

      if (nrow(tmpOut) > 0) {
        if (runWithNormal) {
          tryCatch({
          boxplot(tmpOut[, 2], tmpOut[, 4], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, ylim = c(-2, 2), names = c(HLA_A_type1, HLA_A_type2), las = 1, main = paste('Paired t.test p.val = ', signif (PairedTtest$p.value, 3)), ylab = ('logR ratio'))
          }, error = function(e) return(NULL))
        }
        if (!runWithNormal) {
          tryCatch({
          boxplot(tmpOut[, 2], tmpOut[, 4], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, names = c(HLA_A_type1, HLA_A_type2), las = 1, main = paste('Paired t.test p.val = ', signif (PairedTtest$p.value, 3)), ylab = ('Coverage'))
          }, error = function(e) return(NULL))
        }

        tryCatch({
        beeswarm(tmpOut[, 2], col = c('#de2d2699'), add = TRUE, corral = 'wrap', method = 'swarm', corralWidth = 0.25, pch = 16, at = 1, cex = 1.75)
        }, error = function(e) return(NULL))
        tryCatch({
        beeswarm(tmpOut[, 4], col = c('#3182bd99'), add = TRUE, corral = 'wrap', method = 'swarm', corralWidth = 0.25, pch = 16, at = 2, cex = 1.75)
        }, error = function(e) return(NULL))
        tryCatch({
        barplot(sort(c(tmpOut[, 2]-tmpOut[, 4])), horiz = TRUE, yaxt = 'n', col = ifelse(sort(c(tmpOut[, 2]-tmpOut[, 4]))>0, '#de2d26', '#3182bd'), border = FALSE, main = 'Paired Differences in logR')
        }, error = function(e) return(NULL))
      }

      # t-test of mismatch sites without counting the same read twice
      if (nrow(tmpOut_unique) > 0) {
        if (runWithNormal) {
          tryCatch({
          boxplot(tmpOut_unique[, 2], tmpOut_unique[, 4], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, ylim = c(-10, 10), names = c(HLA_A_type1, HLA_A_type2), las = 1, main = paste('Paired t.test p.val = ', signif (PairedTtest_unique$p.value, 3)), ylab = ('"logR ratio"'))
          }, error = function(e) return(NULL))
        }
        if (!runWithNormal) {
          tryCatch({
          boxplot(tmpOut_unique[, 2], tmpOut_unique[, 4], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, names = c(HLA_A_type1, HLA_A_type2), las = 1, main = paste('Paired t.test p.val = ', signif (PairedTtest_unique$p.value, 3)), ylab = ('"Coverage"'))
          }, error = function(e) return(NULL))
        }
        tryCatch({
        beeswarm(tmpOut_unique[, 2], col = c('#de2d2699'), add = TRUE, corral = 'wrap', method = 'swarm', corralWidth = 0.25, pch = 16, at = 1, cex = 1.75)
        }, error = function(e) return(NULL))
        tryCatch({
        beeswarm(tmpOut_unique[, 4], col = c('#3182bd99'), add = TRUE, corral = 'wrap', method = 'swarm', corralWidth = 0.25, pch = 16, at = 2, cex = 1.75)
        }, error = function(e) return(NULL))
        tryCatch({
        barplot(sort(c(tmpOut_unique[, 2]-tmpOut_unique[, 4])), horiz = TRUE, yaxt = 'n', col = ifelse(sort(c(tmpOut_unique[, 2]-tmpOut_unique[, 4]))>0, '#de2d26', '#3182bd'), border = FALSE, main = 'Paired Differences in logR')
        }, error = function(e) return(NULL))
      }

      # compare tumor / normal coverage site-wise
      # probably not a necessary plot
      if (!any(c(length(HLA_A_type1tumorCov_mismatch_unique),
          length(HLA_A_type2tumorCov_mismatch_unique)) == 0)) {
        tryCatch({
        par(mfrow = c(2, 2))
        tmpOut2 <- cbind(missMatchseq1, MultFactor*HLA_A_type1tumorCov[as.character(missMatchseq1)], HLA_A_type1normalCov[as.character(missMatchseq1)], missMatchseq2, MultFactor*HLA_A_type2tumorCov[as.character(missMatchseq2)], HLA_A_type2normalCov[as.character(missMatchseq2)])
        tmpOut_unique2 <- cbind(missMatchseq1, MultFactor*HLA_A_type1tumorCov_mismatch_unique[as.character(missMatchseq1)], HLA_A_type1normalCov_mismatch_unique[as.character(missMatchseq1)], missMatchseq2, MultFactor*HLA_A_type2tumorCov_mismatch_unique[as.character(missMatchseq2)], HLA_A_type2normalCov_mismatch_unique[as.character(missMatchseq2)])
        dup1_unique <- unique(tmpOut_unique2[duplicated(tmpOut_unique2[, 1]), 1])
        dup2_unique <- unique(tmpOut_unique2[duplicated(tmpOut_unique2[, 4]), 4])

        for (duplicationIn1 in dup1_unique) {
          tmpOut2[tmpOut2[, 1] == duplicationIn1, 5] <-
            median(tmpOut2[tmpOut2[, 1] == duplicationIn1, 5])
          tmpOut2[tmpOut2[, 1] == duplicationIn1, 6] <-
            median(tmpOut2[tmpOut2[, 1] == duplicationIn1, 6])

          tmpOut_unique2[tmpOut_unique2[, 1] == duplicationIn1, 5] <-
            median(tmpOut_unique2[tmpOut_unique2[, 1] == duplicationIn1, 5])
          tmpOut_unique2[tmpOut_unique2[, 1] == duplicationIn1, 6] <-
            median(tmpOut_unique2[tmpOut_unique2[, 1] == duplicationIn1, 6])
        }

        for (duplicationIn2 in dup2_unique) {
          tmpOut2[tmpOut2[, 4] == duplicationIn1, 2] <-
            median(tmpOut2[tmpOut2[, 4] == duplicationIn1, 2])
          tmpOut2[tmpOut2[, 4] == duplicationIn1, 3] <-
            median(tmpOut2[tmpOut2[, 4] == duplicationIn1, 3])

          tmpOut_unique2[tmpOut_unique2[, 4] == duplicationIn2, 2] <-
            median(tmpOut_unique2[tmpOut_unique2[, 4] == duplicationIn2, 2])
          tmpOut_unique2[tmpOut_unique2[, 4] == duplicationIn2, 3] <-
            median(tmpOut_unique2[tmpOut_unique2[, 4] == duplicationIn2, 3])
        }

        tmpOut2 <- tmpOut2[!duplicated(tmpOut2[, 1]), , drop = FALSE]
        tmpOut2 <- tmpOut2[!duplicated(tmpOut2[, 4]), , drop = FALSE]
        tmpOut_unique2 <- tmpOut_unique2[!duplicated(tmpOut_unique2[, 1]), ,
          drop = FALSE]
        tmpOut_unique2 <- tmpOut_unique2[!duplicated(tmpOut_unique2[, 4]), ,
          drop = FALSE]
        }, error = function(e) return(NULL))

        if (nrow(tmpOut2) > 0) {
          tryCatch({
          plot(tmpOut2[, 3], tmpOut2[, 2], xlab = 'normal coverage', ylab = 'tumor coverage', pch = 16, col = '#de2d2699', main = HLA_A_type1)
          abline(a = 0, b = 1)
          }, error = function(e) return(NULL))
          tryCatch({
          plot(tmpOut2[, 6], tmpOut2[, 5], xlab = 'normal coverage',
            ylab = 'tumor coverage', pch = 16, col = '#3182bd99',
            main = HLA_A_type2)
          abline(a = 0, b = 1)
          }, error = function(e) return(NULL))
        }
        if (nrow(tmpOut_unique2) > 0) {
          tryCatch({
          plot(tmpOut_unique2[, 3], tmpOut_unique2[, 2], xlab = 'normal unique coverage', ylab = 'tumor unique coverage', pch = 16, col = '#de2d2699', main = HLA_A_type1)
          abline(a = 0, b = 1)
          }, error = function(e) return(NULL))
          tryCatch({
          plot(tmpOut_unique2[, 6], tmpOut_unique2[, 5], xlab = 'normal unique coverage', ylab = 'tumor unique coverage', pch = 16, col = '#3182bd99', main = HLA_A_type2)
          abline(a = 0, b = 1)
          }, error = function(e) return(NULL))
        }
      }

      if (runWithNormal) {
        ## rolling mean log tumor/normal
        tryCatch({
        par(mfrow = c(1, 1))
        par(mar = c(5, 5, 5, 2))
        plot(c(log2(
          rollmean(HLA_A_type2tumorCov, min(500, length(HLA_A_type2tumorCov)))/
          rollmean(HLA_A_type2normalCov, min(500, length(HLA_A_type2normalCov)))*MultFactor)),
          ylim = c(-2, 2),
          col = '#3182bd', pch = 16, lab = 'HLA genomic position',
          lab = 'Log Ratio',
          ain = c(paste("HLA rolling mean balance", sample)), ex = 0.75)
        points(c(log2(rollmean(HLA_A_type1tumorCov, min(500, length(HLA_A_type1tumorCov)))/rollmean(HLA_A_type1normalCov, min(500, length(HLA_A_type1normalCov)))*MultFactor)), col = '#de2d26', pch = 16, cex = 0.75)
        legend('bottomright', legend = c(HLA_A_type2, HLA_A_type1) ,
               lty = 1, col = c('#3182bd99', '#de2d2699'), bty = 'n', cex = 1, lwd = 3)
        }, error = function(e) return(NULL))

        if (extractNONmismatchReads == TRUE) {
          tryCatch({
          plot(c(1:max(HLA_A_type1normal$V2, HLA_A_type2normal$V2, HLA_A_type2tumor$V2, HLA_A_type1tumor$V2)), ylim = c(-2, 2), col = '#3182bd99', pch = 16, lab = 'HLA genomic position', lab = 'Log Ratio',
            ain = c(paste("HLA raw balance", sample)), ex = 0.75, ype = 'n')
          points(c(HLA_A_type2normal$V2)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], log2(c(HLA_A_type2tumorCov/HLA_A_type2normalCov)*MultFactor)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], col = 'black', bg = '#3182bd', pch = 21, cex = 1)
          points(c(HLA_A_type1normal$V2)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], log2(c(HLA_A_type1tumorCov/HLA_A_type1normalCov)*MultFactor)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], col = 'black', bg = '#de2d26', pch = 21, cex = 1)
          abline(h = 0, lty = 'dashed')

          points(c(HLA_type2normal_nomissmatch$V2), log2(c(HLA_type2tumor_nomissmatchCov/HLA_type2normal_nomissmatchCov)*MultFactor), col = '#3182bd99', pch = 16, cex = 1)
          points(c(HLA_type1normal_nomissmatch$V2), log2(c(HLA_type1tumor_nomissmatchCov/HLA_type1normal_nomissmatchCov)*MultFactor), col = '#de2d2699', pch = 16, cex = 1)
          }, error = function(e) return(NULL))



          # comparing mismatch sites to non-mismatch sites
          tryCatch({
          par(mfrow = c(1, 2))
          par(mar = c(5, 5, 5, 2))

          if (length(HLA_type2tumor_nomissmatchCov) > 1 & length(HLA_A_type2normalCov[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2]) > 1) {
            Ttest <- t.test.NA(log(c((HLA_type2tumor_nomissmatchCov+0.01)/HLA_type2normal_nomissmatchCov)*MultFactor, 2), log(c((HLA_A_type2tumorCov+0.01)/HLA_A_type2normalCov)*MultFactor, 2)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2])
            boxplot(log(c((HLA_type2tumor_nomissmatchCov+0.01)/HLA_type2normal_nomissmatchCov)*MultFactor, 2), log(c((HLA_A_type2tumorCov+0.01)/HLA_A_type2normalCov)*MultFactor, 2)[names(HLA_A_type2normalCov) %in% missMatchPositions$diffSeq2], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, ylim = c(-2, 2), names = c('No mismatch reads', 'Mismatch reads'), las = 1, main = paste(HLA_A_type2, '\n t.test p.val = ', signif (Ttest$p.value, 3)), ylab = ('logR ratio'), xlab = 'coverage from')
          }

          if (length(HLA_type1tumor_nomissmatchCov) > 1 & length(HLA_A_type1normalCov[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1]) > 1) {
            Ttest <- t.test.NA(log(c((HLA_type1tumor_nomissmatchCov+0.01)/HLA_type1normal_nomissmatchCov)*MultFactor, 2), log(c((HLA_A_type1tumorCov+0.01)/HLA_A_type1normalCov)*MultFactor, 2)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1])
            boxplot(log(c((HLA_type1tumor_nomissmatchCov+0.01)/HLA_type1normal_nomissmatchCov)*MultFactor, 2), log(c((HLA_A_type1tumorCov+0.01)/HLA_A_type1normalCov)*MultFactor, 2)[names(HLA_A_type1normalCov) %in% missMatchPositions$diffSeq1], col = c('#de2d2699', '#3182bd99'), boxwex = 0.2, ylim = c(-2, 2), names = c('No mismatch reads', 'Mismatch reads'), las = 1, main = paste(HLA_A_type1, '\n t.test p.val = ', signif(Ttest$p.value, 3)), ylab = ('logR ratio'), xlab = 'coverage from')
          }
          }, error = function(e) return(NULL))
        }
      }
    }
    dev.off()
  }
  ### }}} Plotting step

  ## {{{ Write the output
  fdate <- format(Sys.time(), '%Y%m%d')
  HLAoutLoc <- paste(workDir, '/', patientId, '.',
    minCoverageFilter, '.DNA.HLAlossPrediction_CI.', fdate, fnExt, '.tsv', sep = '')
  write_tsv(HLAoutPut, HLAoutLoc)

  ## Remove redundant rows from output
  # system(glue::glue('grep -v \'^TRUE\' {HLAoutLoc} | sort | uniq | \\
  #     sponge {HLAoutLoc}'))

  if (performIntegerCopyNum) {
    HLABAFsummaryLoc <- paste(workDir, '/', patientId, '.',
      minCoverageFilter, '.DNA.IntegerCPN_CI.', fdate, fnExt, '.tsv', sep = '')
    write_tsv(combinedTable, HLABAFsummaryLoc)
  }
  ## }}}

  ## {{{ Clean up tmp files
  if (cleanUp) {
    cmd <- paste('rm ', workDir, '/', '*tumor*', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*normal*', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*sam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*fastq', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*reads', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*temp*bam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/',
      '*/*chr6region.patient.reference.hlas.csorted.bam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/',
      '*/*chr6region.patient.reference.hlas.csorted.noduplicates.bam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/',
      '*/*chr6region.patient.reference.hlas.bam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*type*[0-9].bam', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm ', workDir, '/', '*/*type*[0-9].bam.bai', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm', workDir, '/', 'mer_counts.jf', sep = '')
    logger(cmd)
    system(cmd)

    cmd <- paste('rm', workDir, '/', 'mer_counts_dumps.fa', sep = '')
    logger(cmd)
    system(cmd)
  }
  ## }}}
}
