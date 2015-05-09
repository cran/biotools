\name{tocher}
\alias{tocher}

\title{ Tocher's Clustering }

\description{
\code{tocher} performs the Tocher (Rao, 1952) optimization clustering from a distance matrix.
}

\usage{
tocher(d, algorithm = c("original", "sequential"))
}

\arguments{
  \item{d}{ an object of class \code{"dist"}.}
  \item{algorithm}{ a character indicating the algorithm to be used for clustering objects.
	It must be one of the two: \code{"original"} (default) or \code{"sequential"}. The latter is 
	the method proposed by Vasconcelos et al. (2007), and sometimes called "modified" Tocher.}
}

\value{
A list of 
  \item{clusters }{a list of length \emph{k} (the number of clusters),
	 containing the labels of the objects in \code{d} for each cluster.}
  \item{class }{a numeric vector indicating the class (the cluster) of each object in \code{d}.}
  \item{criterion }{a numeric vector containing the clustering criteria - the greatest amongst 
	the smallest distances involving each object in \code{d}. If \code{algorithm = "original"},
	this vector contains an unique value, i.e., the same criterion is used for every clustering step.}
  \item{distClust }{a matrix of distances within (diagonal) and between (off-diagonal) clusters.}
  \item{d }{the input object.}
}

\section{Warning }{Clustering a large number of objects (say 300 or more) can be time demanding.}

\references{
Cruz, C.D.; Ferreira, F.M.; Pessoni, L.A. (2011) \emph{Biometria aplicada ao estudo
   da diversidade genetica}. Visconde do Rio Branco: Suprema.

Rao, R.C. (1952) \emph{Advanced statistical methods in biometric research}.
   New York: John Wiley & Sons.

Sharma, J.R. (2006) \emph{Statistical and biometrical techniques in plant breeding}.
   Delhi: New Age International.

Silva, A.R. & Dias, C.T.S. (2013) A cophenetic correlation coefficient for
   Tocher's method. \emph{Pesquisa Agropecuaria Brasileira}, 48:589-596. 

Vasconcelos, E.S.; Cruz, C.D.; Bhering, L.L.; Resende Junior, M.F.R. (2007) Alternative 
   methodology for the cluster analysis. \emph{Pesquisa Agropecuaria Brasileira}, 42:1421-1428.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{dist}}, \code{\link{D2.dist}}, \code{\link{coph.tocher}}, \code{\link{distClust}}, \code{\link{hclust}}
}

\examples{
# example 1
data(garlicdist)
(garlic <- tocher(garlicdist))
garlic$distClust  # cluster distances

# example 2
data(USArrests)
(usa <- tocher(dist(USArrests)))
usa$distClust

# using the sequential algorithm
(usa2 <- tocher(dist(USArrests), algorithm = "sequential"))
usa2$criterion

# example 3
data(eurodist)
(euro <- tocher(eurodist))
euro$distClust

# End (not run)
}
