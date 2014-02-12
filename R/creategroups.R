creategroups <-
function(x, ngroups, sizes, fun = mean,
tol = 1e-2, maxit = 200)
{
   stopifnot(is.numeric(x))
   if (is.null(names(x)))
      stop("'x' must be a named vector!")
   stopifnot(sum(sizes) == length(x))
   stopifnot(ngroups == length(sizes))
   iter <- 0
   repeat {
      g <- sample(rep(1:ngroups, sizes))
      s <- split(x, g)
      label <- split(names(x), g)
      dif <- diff(sapply(s, fun))
      dif <- mean(abs(dif))
      iter <- iter + 1
      if (iter > maxit)
         stop("'maxit' reached!")
      if (dif <= tol)
         break()
   }
   cat("\nCreating homogeneous groups \n")
   cat("------------------------------------------------------------------\n")
   cat("Covariate:",
      deparse(substitute(x)), "\n\nGroups: \n")
   print(label)
   cat("Objective function (equality of):",
      deparse(substitute(fun)), "\n")
   print(sapply(s, fun))
   cat("\nNumber of iterations to convergence:", iter, "\n")
   invisible(s)
}
