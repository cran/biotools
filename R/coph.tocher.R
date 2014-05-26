coph.tocher <-
function(obj)
{
   if (!inherits(obj, "tocher"))
      stop("'obj' must be an object of class 'tocher'!")
   dc <- obj$distClust
   d <- obj$d
   ord <- attr(d, "Labels")
   lab <- unlist(obj$clusters)
   n <- length(lab)
   nc <- length(obj$clusters)
   nopc <- sapply(obj$clusters, length)
   cl <- rep(1:nc, nopc)
   coph <- matrix(0, n, n)
   dimnames(coph) <- list(lab, lab)
   for(j in 1:n) {
      for(i in 1:n) {
         if (i != j) {
            coph[lab[i], lab[j]] <- dc[cl[lab == lab[i]],
               cl[lab == lab[j]]]
         }
      }
   }
   return(as.dist(coph[ord, ord]))
}
