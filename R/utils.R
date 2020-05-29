# This file contains helper functions copied from snpnet


#' @importFrom data.table set as.data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr n
prepareFeatures <- function(pgen, vars, names, stat) {
  buf <- pgenlibr::ReadList(pgen, match(names, vars), meanimpute=F)
  features.add <- as.data.table(buf)
  colnames(features.add) <- names
  for (j in 1:length(names)) {
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stat[["means"]][names[j]])
  }
  features.add
}


computeStats <- function(pfile, ids, configs) {
  keep_f       <- paste0(configs[['gcount.full.prefix']], '.keep')
  gcount_tsv_f <- paste0(configs[['gcount.full.prefix']], '.gcount.tsv')

  dir.create(dirname(configs[['gcount.full.prefix']]), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
      gcount_df <- data.table::fread(gcount_tsv_f)
  } else {      
      # To run plink2 --geno-counts, we write the list of IDs to a file
      data.frame(ID = ids) %>%
      tidyr::separate(ID, into=c('FID', 'IID'), sep='_') %>% 
      data.table::fwrite(keep_f, sep='\t', col.names=F)
  
      # Run plink2 --geno-counts
      system(paste(
          configs[['plink2.path']],
          '--threads', configs[['nCores']],
          '--memory', configs[['mem']],
          '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
          '--keep', keep_f,
          '--out', configs[['gcount.full.prefix']],
          '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs',
          sep=' '
      ), intern=F, wait=T)

      # read the gcount file
      gcount_df <-
        data.table::fread(paste0(configs[['gcount.full.prefix']], '.gcount')) %>%
        dplyr::rename(original_ID = ID) %>%
        dplyr::mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
  }
    
  out <- list()
  out[["pnas"]]  <- gcount_df %>% dplyr::select(stats_pNAs) %>% dplyr::pull()
  out[["means"]] <- gcount_df %>% dplyr::select(stats_means) %>% dplyr::pull()
  out[["sds"]]   <- gcount_df %>% dplyr::select(stats_SDs) %>% dplyr::pull()

  for(key in names(out)){
    names(out[[key]]) <- gcount_df %>% dplyr::select(ID) %>% dplyr::pull()
  }    
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
    
  if (configs[['save']]){
      gcount_df %>% data.table::fwrite(gcount_tsv_f, sep='\t')
      saveRDS(out[["excludeSNP"]], file = file.path(dirname(configs[['gcount.full.prefix']]), "excludeSNP.rda"))
  }

  out
}

readBinMat <- function(fhead, configs){
    # This is a helper function to read binary matrix file (from plink2 --variant-score zs bin)
    rows <- data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', fhead, '.vars.zst'), head=F)$V1
    cols <- data.table::fread(paste0(fhead, '.cols'), head=F)$V1
    bin.reader <- file(paste0(fhead, '.bin'), 'rb')
    M = matrix(
        readBin(bin.reader, 'double', n=length(rows)*length(cols), endian = configs[['endian']]),
        nrow=length(rows), ncol=length(cols), byrow = T
    )
    close(bin.reader)
    colnames(M) <- cols
    rownames(M) <- rows
    if (! configs[['save.computeProduct']]) system(paste(
        'rm', paste0(fhead, '.cols'), paste0(fhead, '.vars.zst'), 
        paste0(fhead, '.bin'), sep=' '
    ), intern=F, wait=T)
    M
}

computeProduct <- function(residual, pfile, vars, stats, configs, iter) {
  time.computeProduct.start <- Sys.time()
  snpnetLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)

  gc_res <- gc()
  if(configs[['KKT.verbose']]) print(gc_res)

  snpnetLogger('Start plink2 --variant-score', indent=3, log.time=time.computeProduct.start)    
  dir.create(file.path(configs[['results.dir']], "save"), showWarnings = FALSE, recursive = T)
    
  residual_f <- file.path(configs[['results.dir']], "save", paste0("residuals_iter_", iter, ".tsv"))
    
  # write residuals to a file
  residual_df <- data.frame(residual)
  colnames(residual_df) <- paste0('lambda_idx_', colnames(residual))
  residual_df %>%    
    tibble::rownames_to_column("ID") %>%
    tidyr::separate(ID, into=c('#FID', 'IID'), sep='_') %>% 
    data.table::fwrite(residual_f, sep='\t', col.names=T)
        
  # Run plink2 --geno-counts
    system(paste(
        configs[['plink2.path']], 
        '--threads', configs[['nCores']],
        '--memory', as.integer(configs[['mem']]) - ceiling(sum(as.matrix(gc_res)[,2])),
        '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
        '--read-freq', paste0(configs[['gcount.full.prefix']], '.gcount'),
        '--keep', residual_f,
        '--out', stringr::str_replace_all(residual_f, '.tsv$', ''),
        '--variant-score', residual_f, 'zs', 'bin', 
        sep=' '
    ), intern=F, wait=T)

  prod.full <- readBinMat(stringr::str_replace_all(residual_f, '.tsv$', '.vscore'), configs)
  if (! configs[['save.computeProduct']] ) system(paste(
      'rm', residual_f, stringr::str_replace_all(residual_f, '.tsv$', '.log'), sep=' '
  ), intern=F, wait=T)
    
  snpnetLoggerTimeDiff('End plink2 --variant-score.', time.computeProduct.start, indent=4)
    
  rownames(prod.full) <- vars    
  if (configs[["standardize.variant"]]) {
      for(residual.col in 1:ncol(residual)){
        prod.full[, residual.col] <- apply(prod.full[, residual.col], 2, "/", stats[["sds"]])
      }
  }
  prod.full[stats[["excludeSNP"]], ] <- NA
  snpnetLoggerTimeDiff('End computeProduct().', time.computeProduct.start, indent=3)
  prod.full
}


setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, covariates, family) {
    if (!("mem" %in% names(configs)))
        stop("mem should be provided to guide the memory capacity.")        
    defaults <- list(
        missing.rate = 0.1, 
        MAF.thresh = 0.001, 
        nCores = 1,
        mem = NULL,
        nlams.init = 10,
        nlams.delta = 5,
        num.snps.batch = 1000, 
        vzs=TRUE, # geno.pfile vzs
        increase.size = NULL,        
        standardize.variant = FALSE,
        early.stopping = TRUE,
        stopping.lag = 2,
        nlambda = 100, 
        niter = 10,
        lambda.min.ratio = NULL,
        glmnet.thresh = 1E-7,
        KKT.verbose = FALSE,
        use.glmnetPlus = NULL,
        save = FALSE,
        save.computeProduct = FALSE,
        prevIter = 0, 
        results.dir = NULL,
        meta.dir = 'meta',
        save.dir = 'results',
        verbose = FALSE,
        KKT.check.aggressive.experimental = FALSE,
        gcount.basename.prefix = 'snpnet.train',
        gcount.full.prefix=NULL,
        endian="little",
        metric=NULL,
        plink2.path='plink2',
        zstdcat.path='zstdcat'
    )
    out <- defaults    
    
    # update the defaults with the specified parameters
    for(name in intersect(names(defaults), names(configs))){
        out[[name]] <- configs[[name]]
    }
    # store additional params
    out[['genotype.pfile']] <- genotype.pfile
    out[['phenotype.file']] <- phenotype.file
    out[['phenotype']] <- phenotype
    out[['covariates']] <- covariates
    out[['family']] <- family
    
    # update settings
    out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
    if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2
    out[['use.glmnetPlus']] <- checkGlmnetPlus(out[['use.glmnetPlus']], family)
        
    if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
    
    # We will write some intermediate files to meta.dir and save.dir.
    # those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
    if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
    dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
    dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
    if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
        out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
    )
    
    out
}
                                 
## logger functions

snpnetLogger <- function(message, log.time = NULL, indent=0, funcname='snpnet'){
    if (is.null(log.time)) log.time <- Sys.time()
    cat('[', as.character(log.time), ' ', funcname, '] ', rep(' ', indent * 2), message, '\n', sep='')
}

timeDiff <- function(start.time, end.time = NULL) {
    if (is.null(end.time)) end.time <- Sys.time()    
    paste(round(end.time-start.time, 4), units(end.time-start.time))
}

snpnetLoggerTimeDiff <- function(message, start.time, end.time = NULL, indent=0){
    if (is.null(end.time)) end.time <- Sys.time()
    snpnetLogger(paste(message, "Time elapsed:", timeDiff(start.time, end.time), sep=' '), log.time=end.time, indent=indent)
}

