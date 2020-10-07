#' @export
solve_aligned <- function(X, y_list, status_list, lambda_1, lambda_2, p.fac = NULL, 
    B0 = NULL)
    {
    K <- length(y_list)
    N <- nrow(X)
    p <- ncol(X)
    order_list <- list()
    status <- matrix(nrow = N, ncol = K)
    rankmin <- matrix(0L, nrow = N, ncol = K)
    rankmax <- matrix(0L, nrow = N, ncol = K)
    for (k in 1:K)
    {
        y <- y_list[[k]]
        o <- order(y)
        y <- y[o]
        order_list[[k]] <- order(o) - 1L
        # status[,k] = status_list[[k]][o]/N # use this vector to control the weight
        status[, k] <- status_list[[k]][o]/sum(status_list[[k]])
        rankmin[, k] <- rank(y, ties.method = "min") - 1L
        rankmax[, k] <- rank(y, ties.method = "min") - 1L
    }
    if (is.null(p.fac))
    {
        p.fac <- rep(1, p)
    }
    if (is.null(B0))
    {
        B0 <- matrix(0, p, K)
    } else
    {
        if (ncol(B0) < K)
        {
            stop("must have at least K columns")
        }
        if (ncol(B0) > K)
        {
            msg <- paste("The number of columns of B0 is greater than K,", "B0[,(K+1):] will be used for proximal step.")
            print(msg)
        }
    }
    
    fit_aligned(X, status, rankmin, rankmax, order_list, B0, lambda_1, lambda_2, 
        p.fac, 1, 5000)
}

#' @export
get_residual <- function(X, y_list, status_list, B)
{
    K <- length(y_list)
    N <- nrow(X)
    p <- ncol(X)
    if (nrow(B) != p | ncol(B) != K)
    {
        stop("The dimension of the parameters is incorrect")
    }
    order_list <- list()
    status <- matrix(nrow = N, ncol = K)
    rankmin <- matrix(0L, nrow = N, ncol = K)
    rankmax <- matrix(0L, nrow = N, ncol = K)
    for (k in 1:K)
    {
        y <- y_list[[k]]
        o <- order(y)
        y <- y[o]
        order_list[[k]] <- order(o) - 1L
        # status[,k] = status_list[[k]][o]/N # use this vector to control the weight
        status[, k] <- status_list[[k]][o]/sum(status_list[[k]])
        rankmin[, k] <- rank(y, ties.method = "min") - 1L
        rankmax[, k] <- rank(y, ties.method = "min") - 1L
    }
    compute_residual(X, status, rankmin, rankmax, order_list, B)
    
}

#' @export
get_dual_norm <- function(grad, alpha, tol = 1e-10)
{   
    if(is.null(dim(grad)))
    {
        grad = as.matrix(grad, ncol=1)
    }
    dnorm <- compute_dual_norm(grad, alpha, tol)
    names(dnorm) <- rownames(grad)
    return(dnorm)
}
