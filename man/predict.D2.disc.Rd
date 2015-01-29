\name{predict.D2.disc}
\alias{predict.D2.disc}

\title{
Classify Multivariate Observations by Mahalanobis Distance 
}

\description{
Classify multivariate observations according to \code{D2.disc} discriminant analysis.
}

\usage{
\method{predict}{D2.disc}(object, newdata, \dots)
}

\arguments{
  \item{object}{an object of class "D2.disc".}
  \item{newdata}{\code{data.frame} or \code{matrix} of observations to be classified.}
  \item{\dots}{arguments based from or to other methods.}
}

\details{
This function is a method for the generic function \code{predict()} for class "D2.disc".
It can be invoked by calling \code{predict(x)} for an object x of the appropriate class,
or directly by calling \code{predict.D2.disc(x)} regardless of the class of the object. 
}

\value{
A list of
  \item{class }{a factor containing the classification of each row in \code{newdata}.}
  \item{D2 }{a matrix (\emph{n x k}) containing Mahalanobis distances between each row in 
	\code{newdata} and the center of each class.}
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\seealso{
\code{\link{D2.disc}}
}

\examples{
disc <- D2.disc(iris[, -5], iris[, 5])
first10 <- iris[1:10, -5]
predict(disc, first10)
predict(disc, iris[, -5])$class

# End (not run)
}
