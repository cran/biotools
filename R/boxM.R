boxM <-
function(data, grouping)
{
   if (!inherits(data, c("data.frame", "matrix")))
      stop("'data' must be a numeric data.frame or matrix!")
   if (length(grouping) != nrow(data))
      stop("incompatible dimensions!")
   data <- as.matrix(data)
   grouping <- as.factor(as.character(grouping))
   p <- ncol(data)
   nlev <- nlevels(grouping)
   lev <- levels(grouping)
   dfs <- tapply(grouping, grouping, length) - 1
   mats <- aux <- list()
   for(i in 1:nlev) {
      mats[[i]] <- cov(data[grouping == lev[i], ])
      aux[[i]] <- mats[[i]] * dfs[i]
   }
   names(mats) <- lev
   pooled <- Reduce("+", aux) / sum(dfs)
   logdet <- log(unlist(lapply(mats, det)))
   minus2logM <- sum(dfs) * log(det(pooled)) - sum(logdet * dfs)
   sum1 <- sum(1 / (dfs - 1)) 
   sum2 <- sum(1 / ((dfs - 1)^2)) 
   Co <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
     (nlev - 1))) * (sum1 - (1 / sum(dfs)))
   X2 <- minus2logM * (1 - Co)
   dfchi <- (choose(p, 2) + p) * (nlev - 1)
   pval <- pchisq(X2, dfchi, lower.tail = FALSE)
   out <- list(cov = mats, pooled = pooled, logDet = logdet,
      chi = X2, df = dfchi, p.value = pval)
   class(out) <- "boxM"
   cat("\n          Box's M-test for Homogeneity of Covariance Matrices \n")
   cat("\nChi-square (approx.) : ", X2, " on ", dfchi,
       "degrees of freedom",
       "\np-value :", pval, "\n")
   invisible(out)
}
