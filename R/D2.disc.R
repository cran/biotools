D2.disc <- 
function(data, grouping, pooled.cov = NULL) 
{
    if (!inherits(data, c("data.frame", "matrix"))) 
        stop("'data' must be a numeric data.frame or matrix!")
    if (length(grouping) != nrow(data)) 
        stop("incompatible dimensions!")
    data <- as.matrix(data)
    name.fac <- deparse(substitute(grouping))
    grouping <- as.factor(as.character(grouping))
    n <- nrow(data)
    p <- ncol(data)
    nlev <- nlevels(grouping)
    lev <- levels(grouping)

    if (is.null(pooled.cov)) {
       dfs <- tapply(grouping, grouping, length) - 1
       if (any(dfs < p)) warning("such a few observations for many variables!")
       mats <- aux <- list()
       for (i in 1:nlev) {
          mats[[i]] <- cov(data[grouping == lev[i], ])
          aux[[i]] <- mats[[i]] * dfs[i]
       }
       names(mats) <- lev
       pooled.cov <- Reduce("+", aux)/sum(dfs)
    } else if (!is.matrix(pooled.cov)) {
       stop("'pooled.cov' must be a square matrix!")
    } else if (dim(pooled.cov)[1] != dim(pooled.cov)[2]) {
       stop("'pooled.cov' must be a square matrix!")
    } else if (any(dim(pooled.cov) != p)) {
       stop("'pooled.cov' has incompatible dimensions with 'data'!")
    }

    med <- aggregate(data, list(grouping), mean)
    med <- as.matrix(med[, -1])
    rownames(med) <- lev

    dists <- matrix(NA, n, nlev, dimnames = list(rownames(data), lev))
    for(i in 1:n) {
       for(j in 1:nlev) {
          dists[i, j] <- mahalanobis(data[i, ], med[j, ], pooled.cov)
       }
    }

    id <- function(x) colnames(dists)[which.min(x)]
    pred <- apply(dists, 1, id)
    misclass <- character(n)
    for(i in 1:n) if (grouping[i] != pred[i]) misclass[i] <- "*"
    confusion <- confusionmatrix(grouping, pred)
    out <- list(D2 = data.frame(dists, grouping, pred, misclass),
       means = med, pooled = pooled.cov, confusion.matrix = confusion)
    class(out) <- "D2.disc"
    return(out)
}
