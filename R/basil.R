#' @importFrom data.table ':='
#' @importFrom data.table set as.data.table
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
#' @export
basil = function(genotype.pfile, phe.file, responsid, covs = NULL, 
                 nlambda = 100, lambda.min.ratio = 0.01, 
                 alpha=NULL, p.factor = NULL,configs = NULL,
                num_lambda_per_iter = 10, num_to_add = 200, use.gpu=FALSE)
{
  if(use.gpu){
    if(exists('solve_aligned_gpu')){
      fit.fun = solve_aligned_gpu
    } else {
      stop("To use our GPU implementation you need to install 'proxgpu' pacakge")
    }
  } else {
    print("To use our GPU implementation, install 'proxgpu' pacakge and set use.gpu=TRUE")
    fit.fun = solve_aligned
  }
  ### Get ids specified by psam --------------------------------------
  psamid = data.table::fread(paste0(genotype.pfile, '.psam'),
                             colClasses = list(character=c("IID")), select = c("IID"))
  psamid = psamid$IID
 
  ### Read responses and covariates --------------------------------------
  status = paste0("coxnet_status_f.", responsid, ".0.0")
  responses = paste0("coxnet_y_f.", responsid, ".0.0")
  
  phe = data.table::fread(phe.file, 
                    colClasses = list(character=c("FID"), factor=c("split")), 
                    select = c("FID", "split", status, responses, covs))
  # Do not allow NA in any column 
  phe=phe[complete.cases(phe), ]
  names(phe)[1] = "ID"
  
  ### Filter out responses with too few events --------------------------------------
  id_to_remove = NULL
  for(i in 1:length(status)){
      s = status[i]
      num_event = sum(phe %>% filter(split == "train") %>% select(all_of(s)))
      cat(paste("Code:", responsid[i]))
      cat(paste("Number of events in validation set", sum(phe %>% filter(split == "val") %>% select(all_of(s)))))
      cat(paste("Number of events in training set", num_event))
      cat("\n")
      if(num_event <100){
          id_to_remove = c(id_to_remove,  responsid[i])        
      }
  }
  
  status_to_remove = paste0("coxnet_status_f.", id_to_remove, ".0.0")
  response_to_remove = paste0("coxnet_y_f.", id_to_remove, ".0.0")
  phe = select(phe, -all_of(c(status_to_remove, response_to_remove)))
  
  status=status[!(responsid %in% id_to_remove)]
  responses = responses[!(responsid %in% id_to_remove)]
  responsid = responsid[!(responsid %in% id_to_remove)]
  
  K = length(responsid) # Number of responses
  if(is.null(alpha)){
    alpha = sqrt(K) # Here alpha is the ratio of lambda_2 and lambda_1
  }
  
  ### Split the data according to the split column ---------------------------------
  phe_train = as.data.table(phe %>% filter(split=='train'))
  phe_val = as.data.table(phe %>% filter(split=='val'))
  
  rm(phe)
  
  ### Initialize train and validation C-index --------------------------------------------
  Ctrain = matrix(-1, nrow=K, ncol=nlambda)
  Cval = matrix(-1, nrow=K, ncol=nlambda)
  
  
  ### Read genotype files, copied from snpnet --------------------------------------------------
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), 
                        VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst'))
  pgen_train = pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(phe_train$ID, psamid))
  pgen_val = pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(phe_val$ID, psamid))


  pgenlibr::ClosePvar(pvar)    

  stats <- computeStats(genotype.pfile, paste(phe_train$ID, phe_train$ID, sep="_"), configs = configs)
  
  ### Fit an unpenalized model ------------------------------------------------------
  if(length(covs) < 1){
    stop("The version without covariates will be implemented later")
  }
  X = as.matrix(select(phe_train, all_of(covs)))
  y_list = list()
  status_list = list()
  for(i in 1:length(responsid)){
    y_list[[i]] = phe_train[[responses[i]]]
    status_list[[i]] = phe_train[[status[i]]]
  }

  result = fit.fun(X,y_list, status_list, c(0.0), c(0.0))
  residuals =  result[['residual']][[1]]
  result = result[['result']]

  ### Compute CIndex ----------------------------------
  X_val = as.matrix(select(phe_val, all_of(covs)))
  for(i in 1:K){
    beta = result[[1]][, i]
    Ctrain[i,1] = cindex::CIndex(X %*% beta, y_list[[i]], status_list[[i]])
    Cval[i,1] = cindex::CIndex(X_val %*% beta, phe_val[[responses[i]]], phe_val[[status[i]]])
  }

  ### Compute residuals and gradient-------------------------------
  residuals = matrix(residuals,nrow = length(phe_train$ID), ncol = K, dimnames = list(paste(phe_train$ID, phe_train$ID, sep='_'), 
                                                                                               paste0("lambda_0_k", 1:K)))

  gradient = computeProduct(residuals, genotype.pfile, vars, stats, configs, iter=0)
  gradient = gradient[-which(rownames(gradient) %in% stats$excludeSNP), ]

  ### Get the dual_norm of the gradient ---------------------------
  score = get_dual_norm(gradient, alpha)
  
  ### Get lambda sequences --------------------------------------------------------
  lambda_max = max(score)
  lambda_min = lambda_max * lambda.min.ratio
  lambda_seq = exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = nlambda))
  # lambda_1 is lamdba_seq, lambda_2 is lambda_seq * alpha
  
  # The first lambda solution is already obtained
  max_valid_index = 1
  prev_valid_index = 0

  # Use validation C-index to determine early stop
  max_cindex = Cval[,1]
  early_stop = rep(FALSE, K)
  names(early_stop) = responsid
  current_response = responsid
  out = list()
  out[[1]] = result[[1]]
  features.to.discard = NULL
  
  iter = 1
  ever.active = covs
  print(ever.active)
  current_B = result[[1]]
  num_not_penalized = length(covs)
  
  ### Start BASIL -----------------------------------------------------------------
  while(max_valid_index < nlambda){

    prev_valid_index = max_valid_index
    print(paste("current maximum valid index is:",max_valid_index ))
    print("Current validation C-Indices are:")
    print(Cval[, 1:max_valid_index])

    if(length(features.to.discard) > 0){
        phe_train[, (features.to.discard) := NULL]
        phe_val[, (features.to.discard) := NULL]
        current_B = current_B[!covs %in% features.to.discard, ,drop=F]
        covs = covs[!covs %in% features.to.discard] 
    }
    
    which.in.model <- which(names(score) %in% covs)
    score[which.in.model] <- NA
    sorted.score <- sort(score, decreasing = T, na.last = NA)

    # if(num_to_add > (16000 - length(covs))){
    #   warning("GPU memory limit will be reached, reduce number of variabels to add")
    #   num_to_add = max(16000 - length(covs),0)
    #   if(num_to_add == 0){
    #     stop("GPU memory limit reached.")
    #   }
    # }
    features.to.add <- names(sorted.score)[1:min(num_to_add, length(sorted.score))]
    covs = c(covs, features.to.add)
    B_init = rbind(current_B, matrix(0.0, nrow=length(features.to.add), ncol=K))
    
    tmp.features.add <- prepareFeatures(pgen_train, vars, features.to.add, stats)
    phe_train[, colnames(tmp.features.add) := tmp.features.add]
    
    tmp.features.add <- prepareFeatures(pgen_val, vars, features.to.add, stats)
    phe_val[, colnames(tmp.features.add) := tmp.features.add]
    
    rm(tmp.features.add)
    
    # Not fit a regularized Cox model for the next few lambdas
    lambda_seq_local = lambda_seq[(max_valid_index + 1):min(max_valid_index + num_lambda_per_iter, length(lambda_seq))]
    # Need better ways to set p.fac
    p.fac = rep(1, nrow(B_init))
    p.fac[1:num_not_penalized] = 0.0
    print(paste("Number of variables to be fitted is:",nrow(B_init)))

    X = as.matrix(select(phe_train, all_of(covs)))

    result = fit.fun(X,y_list, status_list, lambda_seq_local, lambda_seq_local*alpha, p.fac=p.fac, B0=B_init)
    
    residual_all = result[['residual']]
    result = result[['result']]

    residual_all = do.call(cbind, residual_all)
    residual_all = matrix(residual_all,nrow = length(phe_train$ID), ncol = K*num_lambda_per_iter, 
                          dimnames = list(paste(phe_train$ID, phe_train$ID, sep='_'), paste0("lambda_0_k", 1:(K*num_lambda_per_iter))))
    
    gradient = computeProduct(residual_all, genotype.pfile, vars, stats, configs, iter=iter)
    gradient = gradient[-which(rownames(gradient) %in% stats$excludeSNP), ,drop=F]
    
    dnorm_list = list()
    for(i in 1:length(result)){
        start = (i-1)*K+1
        end = i*K
        grad_local = gradient[,start:end, drop=F]
        dnorm_list[[i]] = get_dual_norm(grad_local, alpha)
    }
    

    max_score = sapply(dnorm_list, function(x){max(x[!names(x) %in% covs], na.rm=NA)})
    print("current lambdas are:")
    print(lambda_seq_local)
    print("current Maximum Scores are:")
    print(max_score)
    # if all failed (first one fails)
    if(max_score[1] > lambda_seq_local[1]){
        features.to.discard = NULL
        current_B = result[[1]]
        score = dnorm_list[[1]]
    } else {
        local_valid = which.min(c(max_score <= lambda_seq_local, FALSE)) - 1 # number of valid this iteration
        
        X_val = as.matrix(select(phe_val, all_of(covs)))
        for(j in 1:local_valid){
            out[[max_valid_index+j]] = result[[j]]
            for(i in 1:K){
              beta = result[[j]][, i]
              ind = match(current_response[i], responsid)
              Ctrain[ind,max_valid_index+j] = cindex::CIndex(X %*% beta, y_list[[i]], status_list[[i]])
              Cval[ind,max_valid_index+j] = cindex::CIndex(X_val %*% beta, phe_val[[responses[i]]], phe_val[[status[i]]])
            }
        }
        # Save temp result to files
        save_list = list(Ctrain = Ctrain, Cval = Cval,  beta=out, covs=covs)
        save(save_list, 
                file=file.path(configs[['save.dir']], paste0("saveresult", iter, ".RData")))


        max_Cval_this_iter = apply(Cval[,(max_valid_index + 1):(max_valid_index+local_valid), drop=F], 1, max)
        last_Cval_this_iter = Cval[,(max_valid_index+local_valid)]
        print(max_Cval_this_iter)
        print(last_Cval_this_iter)
        max_cindex = pmax(max_cindex, max_Cval_this_iter)
        #early_stop = early_stop | (signif(last_Cval_this_iter,2) < signif(max_cindex,2))
        early_stop = early_stop | (last_Cval_this_iter < max_cindex - 0.01)

        if(all(early_stop)){
          print("Early stop reached")
          break
        }

        # remove responses whose val Cindex decreases
        if(sum(1-early_stop) < K){
          print("let's remove some response!")
          keep_ind = !(current_response %in% names(which(early_stop)))
          gradient = gradient[,((local_valid - 1)*K+1):(local_valid*K), drop=F]
          gradient = gradient[, keep_ind, drop=F]
          print(keep_ind)
          print(early_stop)
          responses = responses[keep_ind]
          status = status[keep_ind]
          y_list = y_list[keep_ind]
          status_list =  status_list[keep_ind]
          current_B = result[[local_valid]][,keep_ind, drop=F]
          K = sum(1-early_stop)
          alpha = sqrt(K)
          current_response = responsid[!early_stop]
          score = get_dual_norm(gradient, alpha)
          # recompute ever active set, focusing only on the remaning response
          # the ever active set is computed only from valid results in this iteration
          ever.active = covs[unique(unlist(lapply(result, function(x){ which(apply(abs(x), 1, function(y){sum(y[keep_ind])!=0}))})))]
          print(ever.active[1:20])

        } else {
          score =  dnorm_list[[local_valid]]
          current_B = result[[local_valid]]
          new.active = lapply(result, function(x){ which(apply(abs(x), 1, function(y){sum(y)!=0}))})
          ever.active <- union(ever.active, covs[unique(unlist(new.active))])
        }
        print(paste("current number of reponses", K))
        
        max_valid_index = max_valid_index + local_valid
        features.to.discard = setdiff(covs, ever.active)
        print(paste("Number of features discarded in this iteration is", length(features.to.discard)))
        print(paste("Number of ever active variables is", length(ever.active)))
    }
    iter = iter + 1

  }
  return(save_list)
}