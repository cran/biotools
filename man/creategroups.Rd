\name{creategroups}
\alias{creategroups}

\title{ Creating Homogeneous Groups }

\description{
A function to create homogeneous groups of named objects according
to an objective function evaluated at a covariate. It can be useful
to design experiments which contain a fixed covariate factor.
}

\usage{
creategroups(x, ngroups, sizes, fun = mean, tol = 0.01, maxit = 200)
}

\arguments{
  \item{x}{ a numeric vector of a covariate at which to evaluate the objective function.}
  \item{ngroups}{ the number of groups to create.}
  \item{sizes}{ a numeric vector of length equal to \code{ngroups} containing the group sizes. }
  \item{fun}{ the objective function, i.e., to create groups with similar \code{fun}; default is \code{mean}. }
  \item{tol}{ the tolerance level to define the groups as homogenenous; see details.}
  \item{maxit}{ the maximum number of iterations; default is 200. }
}

\details{
\code{creategroups} uses a \code{tol} value to evaluate the following statistic:
\eqn{h = \sum_{j}^{ngroups} abs( t_{j+1} - t_j ) / ngroups}, where \eqn{t_j = fun(group_j)}.
If \eqn{h \leq tol}, the groups are considered homogeneous.
}

\value{
The value returned is a list of named vectors containing the values
for the groups. Also, the result of the objective function evaluated
at the covariate \code{x} is always printed, as well as the number of iterations
required to achieve convergence.
}

\author{
Anderson Rodrigo da Silva <anderson.agro@hotmail.com>
}

\examples{
x <- rnorm(10, 1, 0.5)
names(x) <- letters[1:10]
creategroups(x, ngroups = 2, sizes = c(5, 5))
creategroups(x, ngroups = 3, sizes = c(3, 4, 3), tol = 0.05)

# End (not run)
}
