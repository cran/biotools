coph.tocher <-
function(obj)
{
   if (!inherits(obj, "tocher"))
      stop("'obj' must be an object of class 'dist'!")
   dc <- attr(obj, "distClust")
   d <- attr(obj, "d")
   ord <- attr(d, "Labels")
   lab <- unlist(obj)
   n <- length(lab)
   nc <- length(obj)
   nopc <- sapply(obj, length)
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
