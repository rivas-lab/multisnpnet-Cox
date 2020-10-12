#' @importFrom data.table ':='
#' @importFrom data.table set as.data.table
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
#' @export
basil_base <- function(genotype.pfile, phe.file, responsid, covs, nlambda, lambda.min.ratio, 
    alpha, p.factor, configs, num_lambda_per_iter, num_to_add, max_num_to_add, fit.fun)
    {
    time.start <- Sys.time()
    responsid <- as.character(responsid)
    ### Get ids specified by psam --------------------------------------
    psamid <- data.table::fread(paste0(genotype.pfile, ".psam"), colClasses = list(character = c("IID")), 
        select = c("IID"))
    psamid <- psamid$IID
    
    ### Read responses and covariates --------------------------------------
    time.readphe.start <- Sys.time()
    snpnetLogger("Start reading phenotype file", indent = 2, log.time = time.readphe.start)
    status <- paste0("coxnet_status_f.", responsid, ".0.0")
    responses <- paste0("coxnet_y_f.", responsid, ".0.0")
    
    phe <- data.table::fread(phe.file, colClasses = list(character = c("FID"), factor = c("split")), 
        select = c("FID", "split", status, responses, covs))
    # Do not allow NA in any column
    phe <- phe[complete.cases(phe), ]
    names(phe)[1] <- "ID"
    
    snpnetLoggerTimeDiff("End reading phenotype file.", time.readphe.start, indent = 3)
    
    ### Filter out responses with too few events --------------------------------------
    id_to_remove <- NULL
    for (i in 1:length(status))
    {
        s <- status[i]
        num_event <- sum(phe %>% filter(split == "val") %>% select(all_of(s)))
        printf("Code: %s, number of events in validation set: %-7d", responsid[i], 
            num_event)
        if (num_event < 30)
        {
            id_to_remove <- c(id_to_remove, responsid[i])
            printf("Too few events, removed!")
        }
        printf("\n")
    }

    if(length(id_to_remove) == length(responsid))
    {
        warning("None of the response has sufficient events to proceed")
        return(0)
    }
    
    status_to_remove <- paste0("coxnet_status_f.", id_to_remove, ".0.0")
    response_to_remove <- paste0("coxnet_y_f.", id_to_remove, ".0.0")
    if(!is.null(id_to_remove))
    {
        phe <- select(phe, -all_of(c(status_to_remove, response_to_remove)))
    }
    
    status <- status[!(responsid %in% id_to_remove)]
    responses <- responses[!(responsid %in% id_to_remove)]  # bad name, this is the name of the y column in phe
    responsid <- responsid[!(responsid %in% id_to_remove)]
    names(status) <- responsid
    names(responses) <- responsid
    
    K <- length(responsid)  # Number of responses, this might change
    if (is.null(alpha))
    {
        alpha <- sqrt(K)  # Here alpha is the ratio of lambda_2 and lambda_1
    }
    K0 <- K  # Number of initial response, this won't change
    
    ### Split the data according to the split column ---------------------------------
    phe_test <- phe %>% filter(split == "test")
    phe <- phe %>% filter(split %in% c("train", "val"))
    sprintf("Size of test set is %d.", nrow(phe_test))
    
    sigma <- numeric(length(covs))
    names(sigma) <- covs
    for (cov in covs)
    {
        mu <- mean(phe[[cov]])
        # 0.6 is the sd of many of the snps 
        sigma[cov] <- sd(phe[[cov]]) / 0.6
        phe[[cov]] <- (phe[[cov]] - mu)/sigma[cov]
        phe_test[[cov]] <- (phe_test[[cov]] - mu)/sigma[cov]
    }
    phe_train <- as.data.table(phe %>% filter(split == "train"))
    phe_val <- as.data.table(phe %>% filter(split == "val"))
    
    rm(phe)
    
    ### Initialize train and validation C-index
    ### --------------------------------------------
    Ctrain <- matrix(-1, nrow = K, ncol = nlambda)
    Cval <- matrix(-1, nrow = K, ncol = nlambda)
    Ctest <- matrix(-1, nrow = K, ncol = nlambda)
    rownames(Ctrain) <- responsid
    rownames(Cval) <- responsid
    rownames(Ctest) <- responsid
    
    
    ### Read genotype files, copied from snpnet
    ### --------------------------------------------------
    time.computestats.start <- Sys.time()
    snpnetLogger("Start loading and computing genotype statistics", indent = 2, log.time = time.computestats.start)
    
    vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd = paste0(configs[["zstdcat.path"]], 
        " ", paste0(genotype.pfile, ".pvar.zst"))), CHROM = "#CHROM"), VAR_ID = paste(ID, 
        ALT, sep = "_"))$VAR_ID
    pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, ".pvar.zst"))
    pgen_train <- pgenlibr::NewPgen(paste0(genotype.pfile, ".pgen"), pvar = pvar, 
        sample_subset = match(phe_train$ID, psamid))
    pgen_val <- pgenlibr::NewPgen(paste0(genotype.pfile, ".pgen"), pvar = pvar, sample_subset = match(phe_val$ID, 
        psamid))

    pgen_test <- pgenlibr::NewPgen(paste0(genotype.pfile, ".pgen"), pvar = pvar, sample_subset = match(phe_test$ID, psamid))
    
    
    pgenlibr::ClosePvar(pvar)
    
    stats <- computeStats(genotype.pfile, paste(phe_train$ID, phe_train$ID, sep = "_"), 
        configs = configs)
    
    snpnetLoggerTimeDiff("End loading and computing genotype statistics.", time.computestats.start, 
        indent = 3)
    
    cov.factor <- rep(0, length(covs))
    names(cov.factor) <- covs
    if (is.null(p.factor))
    {
        p.factor <- rep(1, length(vars))
        names(p.factor) <- vars
    } else
    {
        if (!(all(vars %in% names(p.factor))))
        {
            stop("p.factor must be NULL or provide factor for all variants")
        }
        if (any(p.factor <= 0))
        {
            stop("We only support stricly positive penalty factor.")
        }
    }
    p.factor <- c(cov.factor, p.factor)
    
    ### Fit an unpenalized model ------------------------------------------------------
    if (length(covs) < 1)
    {
        stop("The version without covariates will be implemented later")
    }
    iter <- 0
    # time.basilload.start <- Sys.time() snpnetLogger(sprintf('Start loading training
    # data for basil iteration %d.', iter), indent=2, log.time=time.basilload.start)
    X <- as.matrix(select(phe_train, all_of(covs)))
    y_list <- list()
    status_list <- list()
    for (i in 1:length(responsid))
    {
        y_list[[responsid[i]]] <- phe_train[[responses[i]]]
        status_list[[responsid[i]]] <- phe_train[[status[i]]]
    }
    # snpnetLoggerTimeDiff(sprintf('End loading training data for basil iteration
    # %d.', iter), time.basilload.start, indent=3)
    
    
    time.basilfit.start <- Sys.time()
    snpnetLogger(sprintf("Start fitting step for basil iteration %d.", iter), indent = 2, 
        log.time = time.basilfit.start)
    result <- fit.fun(X, y_list, status_list, c(0), c(0))
    snpnetLoggerTimeDiff(sprintf("End fitting step for basil iteration %d.", iter), 
        time.basilfit.start, indent = 3)
    
    residuals <- result[["residual"]][[1]]
    result <- result[["result"]]
    
    ### Compute CIndex ----------------------------------
    time.basilmetric.start <- Sys.time()
    snpnetLogger(sprintf("Start metric evaluations for basil iteration %d.", iter), 
        indent = 2, log.time = time.basilmetric.start)
    X_val <- as.matrix(select(phe_val, all_of(covs)))
    X_test <- as.matrix(select(phe_test, all_of(covs)))
    pred_train <- X %*% result[[1]]
    pred_val <- X_val %*% result[[1]]
    pred_test <- X_test %*% result[[1]]
    for (i in 1:K)
    {
        Ctrain[i, 1] <- cindex::CIndex(pred_train[, i], y_list[[i]], status_list[[i]])
        Cval[i, 1] <- cindex::CIndex(pred_val[, i], phe_val[[responses[i]]], phe_val[[status[i]]])
        Ctest[i, 1] <- cindex::CIndex(pred_test[, i], phe_test[[responses[i]]], phe_test[[status[i]]])
    }
    snpnetLoggerTimeDiff(sprintf("End metric evaluations for basil iteration %d.", 
        iter), time.basilmetric.start, indent = 3)
    
    ### Compute residuals and gradient-------------------------------
    residuals <- matrix(residuals, nrow = length(phe_train$ID), ncol = K, dimnames = list(paste(phe_train$ID, 
        phe_train$ID, sep = "_"), paste0("lambda_0_k", 1:K)))
    
    gradient <- computeProduct(residuals, genotype.pfile, vars, stats, configs, iter = 0)
    gradient <- gradient[-which(rownames(gradient) %in% stats$excludeSNP), ]
    
    ### Get the dual_norm of the gradient ---------------------------
    time.dnorm.start <- Sys.time()
    snpnetLogger(sprintf("Computing dual-norm for basil iteration %d.", iter), indent = 2, 
        log.time = time.dnorm.start)
    score <- get_dual_norm(gradient, alpha)
    score <- score/p.factor[names(score)]
    snpnetLoggerTimeDiff(sprintf("Finished computing dual-norm for basil iteration %d.", 
        iter), time.dnorm.start, indent = 3)
    
    ### Get lambda sequences --------------------------------------------------------
    lambda_max <- max(score)
    lambda_min <- lambda_max * lambda.min.ratio
    lambda_seq <- exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = nlambda))
    # lambda_1 is lamdba_seq, lambda_2 is lambda_seq * alpha
    
    # The first lambda solution is already obtained
    max_valid_index <- 1
    
    # Use validation C-index to determine early stop
    max_cindex <- Cval[, 1]
    early_stop <- rep(FALSE, K)
    names(early_stop) <- responsid
    current_B <- result[[1]]
    rownames(current_B) <- covs
    colnames(current_B) <- responsid
    out <- list()
    out[[1]] <- current_B
    features.to.discard <- NULL
    
    iter <- 1
    ever.active <- covs
    best_lam_ind <- rep(1, K)  # keep track of iter at which each response achieves best validation metric
    names(best_lam_ind) <- responsid

    current_response <- responsid
    num_not_penalized <- length(covs)
    
    KKT_failure_count <- 0
    ### Start BASIL -----------------------------------------------------------------
    while (max_valid_index < nlambda)
    {
        time.basil.start <- Sys.time()
        snpnetLogger(sprintf("Start basil iteration %d.", iter), indent = 2, log.time = time.basil.start)
        
        prev_valid_index <- max_valid_index
        printf("Current maximum valid index is: %d\n", max_valid_index)
        printf("Current test C-Indices are:\n")
        print(Ctest[, 1:max_valid_index])
        
        if (length(features.to.discard) > 0)
        {
            phe_train[, `:=`((features.to.discard), NULL)]
            phe_val[, `:=`((features.to.discard), NULL)]
            phe_test[, `:=`((features.to.discard), NULL)]
            current_B <- current_B[!covs %in% features.to.discard, , drop = F]
            covs <- covs[!covs %in% features.to.discard]
        }
        
        which.in.model <- which(names(score) %in% covs)
        score[which.in.model] <- NA
        sorted.score <- sort(score, decreasing = T, na.last = NA)
        
        features.to.add <- names(sorted.score)[1:min(num_to_add, length(sorted.score))]
        covs <- c(covs, features.to.add)
        B_init <- rbind(current_B, matrix(0, nrow = length(features.to.add), ncol = ncol(current_B)))
        
        time.preparefeatures.start <- Sys.time()
        snpnetLogger(sprintf("Preparing features for basil iteration %d.", iter), 
            indent = 2, log.time = time.preparefeatures.start)
        tmp.features.add <- prepareFeatures(pgen_train, vars, features.to.add, stats)
        phe_train[, `:=`(colnames(tmp.features.add), tmp.features.add)]
        
        tmp.features.add <- prepareFeatures(pgen_val, vars, features.to.add, stats)
        phe_val[, `:=`(colnames(tmp.features.add), tmp.features.add)]

        tmp.features.add <- prepareFeatures(pgen_test, vars, features.to.add, stats)
        phe_test[, `:=`(colnames(tmp.features.add), tmp.features.add)]
        snpnetLoggerTimeDiff(sprintf("End preparing features for basil iteration %d.", 
            iter), time.preparefeatures.start, indent = 3)
        
        rm(tmp.features.add)
        
        # Not fit a regularized Cox model for the next few lambdas
        lambda_seq_local <- lambda_seq[(max_valid_index + 1):min(max_valid_index + 
            num_lambda_per_iter, length(lambda_seq))]
        
        # p.fac = rep(1, nrow(B_init)) p.fac[1:num_not_penalized] = 0.0
        
        p.fac <- p.factor[covs]
        printf("Number of variables to be fitted is: %d. \n", nrow(B_init))
        
        X <- as.matrix(select(phe_train, all_of(covs)))
        
        time.basilfit.start <- Sys.time()
        snpnetLogger(sprintf("Start fitting step for basil iteration %d.", iter), 
            indent = 2, log.time = time.basilfit.start)
        result <- fit.fun(X, y_list, status_list, lambda_seq_local, lambda_seq_local * 
            alpha, p.fac = p.fac, B0 = B_init)
        snpnetLoggerTimeDiff(sprintf("End fitting step for basil iteration %d.", 
            iter), time.basilfit.start, indent = 3)
        
        residual_all <- result[["residual"]]
        result <- result[["result"]]
        for (l in 1:length(result))
        {
            rownames(result[[l]]) <- covs
            colnames(result[[l]]) <- responsid
        }
        # print(c('Colnames of results:', colnames(result[[1]]))) print(c('Colnames of
        # current_B:', colnames(current_B)))
        
        residual_all <- do.call(cbind, residual_all)
        residual_all <- matrix(residual_all, nrow = length(phe_train$ID), ncol = K * 
            num_lambda_per_iter, dimnames = list(paste(phe_train$ID, phe_train$ID, 
            sep = "_"), paste0("lambda_0_k", 1:(K * num_lambda_per_iter))))
        
        gradient <- computeProduct(residual_all, genotype.pfile, vars, stats, configs, 
            iter = iter)
        gradient <- gradient[-which(rownames(gradient) %in% stats$excludeSNP), , 
            drop = F]
        
        time.dnorm.start <- Sys.time()
        snpnetLogger(sprintf("Computing dual-norm for basil iteration %d.", iter), 
            indent = 2, log.time = time.dnorm.start)
        dnorm_list <- list()
        for (i in 1:length(result))
        {
            start <- (i - 1) * K + 1
            end <- i * K
            grad_local <- gradient[, start:end, drop = F]
            dnorm_list[[i]] <- get_dual_norm(grad_local, alpha)/p.factor[rownames(grad_local)]
        }
        snpnetLoggerTimeDiff(sprintf("End computing dual-norm for basil iteration %d.", 
            iter), time.dnorm.start, indent = 3)
        
        
        max_score <- sapply(dnorm_list, function(x)
        {
            max(x[!names(x) %in% covs], na.rm = NA)
        })
        printf("current lambdas are:\n")
        print(lambda_seq_local)
        printf("current Maximum Scores are:\n")
        print(max_score)
        # if all failed (first one fails)
        if (max_score[1] > lambda_seq_local[1])
        {
            features.to.discard <- NULL
            current_B <- result[[1]]
            score <- dnorm_list[[1]]
            KKT_failure_count <- KKT_failure_count + 1
        } else
        {
            local_valid <- which.min(c(max_score <= lambda_seq_local, FALSE)) - 1  # number of valid this iteration
            
            if (local_valid <= (num_lambda_per_iter/2))
            {
                KKT_failure_count <- KKT_failure_count + 1
            }
            
            time.basilmetric.start <- Sys.time()
            snpnetLogger(sprintf("Start metric evaluations for basil iteration %d.", 
                iter), indent = 2, log.time = time.basilmetric.start)
            X_val <- as.matrix(select(phe_val, all_of(covs)))
            X_test <- as.matrix(select(phe_test, all_of(covs)))
            for (j in 1:local_valid)
            {
                out[[max_valid_index + j]] <- result[[j]]
                for (i in 1:K)
                {
                  beta <- result[[j]][, i]
                  ind <- match(current_response[i], responsid)
                  Ctrain[ind, max_valid_index + j] <- cindex::CIndex(X %*% beta, 
                    y_list[[i]], status_list[[i]])
                  cval_tmp <- cindex::CIndex(X_val %*% beta, phe_val[[responses[i]]], 
                    phe_val[[status[i]]])
                  Cval[ind, max_valid_index + j] <- cval_tmp
                  if (cval_tmp > max_cindex[ind])
                  {
                    max_cindex[ind] <- cval_tmp
                    best_lam_ind[ind] <- max_valid_index + j
                  }
                  Ctest[ind, max_valid_index + j] <- cindex::CIndex(X_test %*% beta, phe_test[[responses[i]]], 
                    phe_test[[status[i]]])




                }
            }
            snpnetLoggerTimeDiff(sprintf("End metric evaluations for basil iteration %d.", 
                iter), time.basilmetric.start, indent = 3)
            # Save temp result to files
            save_list <- list(Ctrain = Ctrain, Cval = Cval, Ctest=Ctest, beta = out)
            save(save_list, file = file.path(configs[["save.dir"]], paste0("saveresult", 
                iter, ".RData")))
            
            last_Cval_this_iter <- Cval[, (max_valid_index + local_valid)]
            # max_Cval_this_iter = apply(Cval[,(max_valid_index +
            # 1):(max_valid_index+local_valid), drop=F], 1, max) early_stop = early_stop |
            # (last_Cval_this_iter < max_cindex)
            early_stop <- early_stop | (last_Cval_this_iter < max_cindex)
            
            # Don't stop too early
            if (all(early_stop) && length(ever.active) > 3000)
            {
                snpnetLoggerTimeDiff("Early stop for all responses reached.", time.start, 
                  indent = 3)
                break
            }
            
 
            score <- dnorm_list[[local_valid]]
            current_B <- result[[local_valid]]
            new.active = NULL
            for (j in 1:local_valid)
            {
                new.active  <- union(new.active, names(which(apply(abs(result[[j]]), 1, function(y){sum(y)!=0}))))
            }
            ever.active <- union(ever.active, new.active)
            
            printf("size of ever active set is %d.\n", length(ever.active))
            print(ever.active[1:20])
            
            max_valid_index <- max_valid_index + local_valid
            features.to.discard <- setdiff(covs, ever.active)
            printf("Number of inactive features discarded in this iteration is %d.\n", 
                length(features.to.discard))
        }
        if (KKT_failure_count >= 2)
        {
            num_to_add <- min(num_to_add * 2, max_num_to_add)
            KKT_failure_count <- 0
        }
        
        snpnetLoggerTimeDiff(sprintf("End basil iteration %d.", iter), time.basil.start, 
            indent = 3)
        iter <- iter + 1
        
    }
    return(save_list)
}


#' @export
basil <- function(genotype.pfile, phe.file, responsid, covs = NULL, nlambda = 100, 
    lambda.min.ratio = 0.01, alpha = NULL, p.factor = NULL, configs = NULL, num_lambda_per_iter = 10, 
    num_to_add = 1500, max_num_to_add = 6000)
    {
    basil_base(genotype.pfile, phe.file, responsid, covs, nlambda, lambda.min.ratio, 
        alpha, p.factor, configs, num_lambda_per_iter, num_to_add, max_num_to_add, 
        solve_aligned)
}
