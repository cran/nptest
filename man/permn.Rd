\name{permn}
\alias{permn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate All Permutations of n Elements
}
\description{
Generates all \eqn{n!} vectors of length \eqn{n} consisting of permutations of the integers 1 to \eqn{n}.
}
\usage{
permn(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of elements.
}
}
\details{
Adapted from the "permutations" function in the \href{https://CRAN.R-project.org/package=e1071}{e1071} R package.
}
\value{
Matrix of dimension \eqn{n} by \eqn{n!} where each column contains a unique permutation vector.
}
\references{
Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., & Leisch, F. (2018). e1071: Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R package version 1.7-0. https://CRAN.R-project.org/package=e1071
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\note{
Used for exact tests in \code{\link{np.cor.test}} and \code{\link{np.reg.test}}.
}

\section{Warning }{
For large \code{n} this function will consume a lot of memory and may even crash R.
}

\examples{
permn(2)
permn(3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
