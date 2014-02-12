\name{tocher}
\alias{tocher}

\title{ Tocher's Clustering }

\description{
\code{tocher} performs the Tocher (Rao, 1952) optimization clustering
from a distance matrix.
}

\usage{
tocher(d)
}

\arguments{
  \item{d}{ an object of class "dist".}
}

\value{
\code{tocher} returns a list of length \emph{k} (the number of clusters) containing
the labels of the objects in \code{d}. Furthermore, a matrix of distances
within (diagonal) and between (off-diagonal) clusters is attached as an attribute.
}

\references{
Cruz, C.D.; Ferreira, F.M.; Pessoni, L.A. (2011) \emph{Biometria aplicada ao estudo
   da diversidade gen\'etica}. Visconde do Rio Branco: Suprema.

Rao, R.C. (1952) \emph{Advanced statistical methods in biometric research}.
   New York: John Wiley & Sons.

Sharma, J.R. (2006) \emph{Statistical and biometrical techniques in plant breeding}.
   Delhi: New Age International.

Silva, A.R. & Dias, C.T.S. (2013) A cophenetic correlation coefficient for
   Tocher's method. \emph{Pesquisa Agropecu\'aria Brasileira}, 48:589-596. 
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{coph.tocher}}, \code{\link{distClust}}, \code{\link{hclust}}
}

\examples{
# example 1
data(garlicdist)
garlic <- tocher(garlicdist)
attr(garlic, "distClust") # cluster distances

# example 2
data(USArrests)
usa <- tocher(dist(USArrests))
attr(usa, "distClust")

# example 3
data(eurodist)
euro <- tocher(eurodist)
attr(euro, "distClust")

# End (not run)
}
