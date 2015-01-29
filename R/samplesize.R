samplesize <- 
function(x, fun, sizes = NULL, lcl = NULL, ucl = NULL, 
	nboot = 200, conf.level = 0.95, nrep = 500, graph = TRUE, ...)
{
    stopifnot(is.numeric(x))
    x <- as.vector(x)
    stopifnot(is.function(fun))
    stopifnot(nboot > 1)
    stopifnot(nrep > 1)
    n <- length(x)

    if (is.null(lcl) || is.null(ucl)) {
       bo <- boot(x, fun, nboot)
       ci = boot.ci(bo, conf.level, type = "perc")
       lcl <- ci$percent[4L]
       ucl <- ci$percent[5L]
    }

    if (is.null(sizes)) sizes <- 2:(n-1)
    aux <- matrix(NA, nrep, length(sizes))
    for(j in 1:length(sizes)) {
       for(i in 1:nrep) {
          aux[i, j] <- fun(sample(x, sizes[j], replace = TRUE))
       }
    }

    f.out <- function(x) length(c(x[x < lcl], x[x > ucl]))
    n.out <- apply(aux, 2, f.out)
    prop <- n.out/nrep
    
    if(graph) {
       plot(aux[1, ] ~ sizes,
          ylim = range(c(aux)),
          ylab = deparse(substitute(fun)),
          xlab = "Sample size", ...)
       for(j in 2:nrep) points(aux[j, ] ~ sizes, ...)
       abline(h = c(lcl, ucl), col = "red")
       devAskNewPage(ask = TRUE)
       plot(n.out ~ sizes,
          type = "b",
          xlab = "Sample size",
          ylab = "N outside the CI", ...)
    }

    out <- list(CI = c(lcl, ucl),
       pointsOut = data.frame(sizes, n.out, prop))
    class(out) <- "samplesize"
    return(out)
}
