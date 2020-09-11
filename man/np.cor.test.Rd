\name{np.cor.test}
\alias{np.cor.test}

\title{
Nonparametric Tests of Correlation Coefficients
}
\description{
Denoting the Pearson product-moment correlation coefficient as
\deqn{\rho = Cov(X, Y) / \sqrt{Var(X) Var(Y)}} 
this function implements permutation tests of \eqn{H_0: \rho = \rho_0} where \eqn{\rho_0} is the user-specified null value. Can also implement tests of partial correlations, semi-partial (or part) correlations, and independence.
}

\usage{
np.cor.test(x, y, z = NULL,
            alternative = c("two.sided", "less", "greater"),
            rho = 0, independent = FALSE, partial = TRUE,
            R = 9999, parallel = FALSE, cl = NULL,
            perm.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\eqn{X} vector (\code{n} by 1).
}
  \item{y}{
\eqn{Y} vector (\code{n} by 1).
}
  \item{z}{
Optional \eqn{Z} matrix (\code{n} by \code{q}). If provided, the partial (or semi-partial if \code{partial = FALSE}) correlation is calculated between \code{x} and \code{y} controlling for \code{z}.
}
  \item{alternative}{
Alternative hypothesis. Must be either "two.sided" (\eqn{H_1: \rho \neq \rho_0}), "less" (\eqn{H_1: \rho < \rho_0}), or "greater" (\eqn{H_1: \rho > \rho_0}).
}
  \item{rho}{
Null hypothesis value \eqn{\rho_0}. Defaults to zero.
}
  \item{independent}{
If \code{FALSE} (default), the null hypothesis is \eqn{H_0: \rho = \rho_0}. Otherwise, the null hythpothesis is that \eqn{X} and \eqn{Y} are independent, i.e., \eqn{H_0: F_{XY}(x,y) = F_X(x) F_Y(y)}.
}
  \item{partial}{
Only applicable if \code{z} is provided. If \code{TRUE} (default), the partial correlation between \code{x} and \code{y} controlling for \code{z} is tested. Otherwise the semi-partial correlation is tested. See Details. 
}
  \item{R}{
Number of resamples for the permutation test (positive integer). 
}
  \item{parallel}{
Logical indicating if the \code{\link{parallel}} package should be used for parallel computing (of the permutation distribution). Defaults to \code{FALSE}, which implements sequential computing.
}
  \item{cl}{
Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}. 
}
  \item{perm.dist}{
Logical indicating if the permutation distribution should be returned.  
}
}
\details{
Default use of this function tests the Pearson correlation between \eqn{X} and \eqn{Y} using the studentized test statistic proposed by DiCiccio and Romano (2017). If \code{independent = TRUE}, the classic (unstudentized) test statistic is used to test the null hypothesis of independence.

If \eqn{Z} is provided, the partial or semi-partial correlation between \eqn{X} and \eqn{Y} controlling for \eqn{Z} is tested. For the semi-partial correlation, the effect of \eqn{Z} is partialled out of \eqn{X}. 
}

\value{
  \item{statistic }{
Test statistic value.
}
  \item{p.value }{
p-value for testing \eqn{H_0: \rho = \rho_0} or \eqn{H_0: F_{XY}(x,y) = F_X(x) F_Y(y)}.
}
  \item{perm.dist }{
Permutation distribution of \code{statistic}. 
}
  \item{alternative }{
Alternative hypothesis.
}
  \item{null.value }{
Null hypothesis value for \eqn{\rho}.
}
  \item{independent }{
Independence test?
}
  \item{R }{
Number of resamples.
}
  \item{exact }{
Exact permutation test? See Note.
}
  \item{estimate }{
Sample estimate of correlation coefficient \eqn{\rho}.
}
}
\references{
DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for correlation and regression coefficients. Journal of the American Statistical Association, 112(519), 1211-1220. doi: 10.1080/01621459.2016.1202117

Helwig, N. E. (2019). Statistical nonparametric mapping: Multivariate permutation tests for location, correlation, and regression problems in neuroimaging. WIREs Computational Statistics, 11(2), e1457. doi: 10.1002/wics.1457

Pitman, E. J. G. (1937b). Significance tests which may be applied to samples from any populations. ii. the correlation coefficient test. Supplement to the Journal of the Royal Statistical Society, 4(2), 225-232. doi: 10.2307/2983647
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\note{
The permutation test will be exact when the requested number of resamples \code{R} is greater than \code{factorial(n)} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{factorial(n)} possible values of the test statistic.

If \code{z = NULL}, the result will be the same as using \code{\link{np.reg.test}} with \code{method = "perm"}.

If \code{z} is supplied and \code{partial = TRUE}, the result will be the same as using \code{\link{np.reg.test}} with \code{method = "KC"} and \code{homosced = FALSE}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.np.cor.test}} S3 plotting method for visualizing the results
}

\examples{
# generate data
rho <- 0.5
val <- c(sqrt(1 + rho), sqrt(1 - rho))
corsqrt <- matrix(c(val[1], -val[2], val), 2, 2) / sqrt(2)
set.seed(1)
n <- 10
z <- cbind(rnorm(n), rnorm(n)) \%*\% corsqrt
x <- z[,1]
y <- z[,2]

# test H0: rho = 0
set.seed(0)
np.cor.test(x, y)

# test H0: X and Y are independent
set.seed(0)
np.cor.test(x, y, independent = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ nonparametric }
\keyword{ robust }
\keyword{ univar }
