\name{pathanalysis}
\alias{pathanalysis}

\title{ Path Analysis, Simple and Under Colinearity }

\description{
Function to perform the simple path analysis and the path analysis
under colinearity (sometimes called \emph{ridge path analysis}). It computes
the direct (diagonal) and indirect (off-diagonal) effects of each explanatory
variable over a response one.
}

\usage{
pathanalysis(corMatrix, resp.col, colinearity = FALSE)
}

\arguments{
  \item{corMatrix}{ a correlation matrix. }
  \item{resp.col}{ an integer value indicating the column in \code{corMatrix} that corresponds to
	the response variable. }
  \item{colinearity}{ logical; if \code{TRUE}, an external interactive display is
	used to pass a value, say k, at which to evaluate the system:
	(\bold{X'X} + \bold{I}k)\bold{B} = \bold{X'Y}, being \bold{X'X} the correlation matrix between explanatory variables,
	\bold{X'Y} the correlation vector between all explanatory variables and the response variable,
	\bold{B} is the vector of path coefficients and \emph{k} is a value between 0 and 1; default is \code{FALSE}, i.e., k = 0. }
}

\value{
A list of
  \item{coef }{ a matrix containing the direct (diagonal) and indirect (off-diagonal)
        effects of each variable.}
  \item{Rsq }{ the coefficient of determination.}
  \item{ResidualEffect }{ the residual effect.}
  \item{VIF }{ a vector containing the variance inflation factors.}
  \item{CN }{ the condition number.}
}

\references{
Carvalho, S.P. (1995) \emph{Metodos alternativos de estimacao de coeficientes de
   trilha e indices de selecao, sob multicolinearidade}. Ph.D. Thesis, 
   Federal University of Vicosa (UFV), Vicosa, MG, Brazil.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\note{
If \code{colinearity} is \code{TRUE}, the outcome is just printed on R console.
}

\examples{
data(peppercorr)
pathanalysis(peppercorr, 6, colinearity = TRUE)

# End (not run)
}
