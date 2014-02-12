\name{boxM}
\alias{boxM}

\title{ Box's M-test }

\description{
It performs the Box's M-test for homogeneity of covariance matrices
obtained from multivariate normal data according to one classification
factor. The test is based on the chi-square approximation.
}

\usage{
boxM(data, grouping)
}

\arguments{
  \item{data}{ a numeric data.frame or matrix containing n observations of p variables;
	it is expected that \code{n > p}. }
  \item{grouping}{ a vector of length n containing the class of each observation;
	it is usualy a factor. }
}

\value{
A list of
  \item{cov }{a list containing the within covariance matrix for each level of \code{grouping}.}
  \item{pooled }{the pooled covariance matrix.}
  \item{logDet }{a vector containing the natural logarithm of each matrix in \code{cov}.}
  \item{chi }{an approximated value of the chi-square distribution.}
  \item{df }{the number of degrees of freedom related to \code{chi}.}
 \item{p.value }{the p-value of the test.}
}

\references{
Morrison, D.F. (1976) \emph{Multivariate Statistical Methods}.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\examples{
data(iris)
boxM(iris[, -5], iris[, 5])

# End (not run)
}
