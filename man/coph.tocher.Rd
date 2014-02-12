\name{coph.tocher}
\alias{coph.tocher}

\title{ Cophenetic Distances for Tocher's Clustering }

\description{
It computes the cophenetic distance matrix for a Tocher's clustering.
}

\usage{
coph.tocher(obj)
}

\arguments{
  \item{obj}{ an object of class \code{\link{tocher}}.}
}


\value{
An object of class \code{dist}.
}

\references{
Silva, A.R. & Dias, C.T.S. (2013) A cophenetic correlation coefficient for
   Tocher's method. \emph{Pesquisa Agropecuaria Brasileira}, 48:589-596. 
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{tocher}}
}

\examples{
# Example 1
data(garlicdist)
garlic <- tocher(garlicdist)
coph <- coph.tocher(garlic)
cor(coph, garlicdist)

# Example 2
usa <- tocher(dist(USArrests))
cophUS <- coph.tocher(usa)
cor(cophUS, dist(USArrests))

# End (not run)
}
