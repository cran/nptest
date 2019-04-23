\name{mcse}
\alias{mcse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Standard Errors for Tests
}
\description{
This function calculates Monte Carlo standard errors for (non-exact) nonparametric tests. The MCSEs can be used to determine (i) the accuracy of a test for a given number of resamples, or (ii) the number of resamples needed to achieve a test with a given accuracy.
}
\usage{
mcse(R, delta, conf.level = 0.95, sig.level = 0.05,
     alternative = c("two.sided", "one.sided"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
Number of resamples (positive integer).
}
  \item{delta}{
Accuracy of the approximation (number between 0 and 1).
}
  \item{conf.level}{
Confidence level for the approximation (number between 0 and 1).
}
  \item{sig.level}{
Significance level of the test (number between 0 and 1).
}
  \item{alternative}{
Alternative hypothesis (two-sided or one-sided).
}
}
\details{
Note: either \code{R} or \code{delta} must be provided.

Let \eqn{F(x)} denote the distribution function for the full permutation distribution, and let \eqn{G(x)} denote the approximation obtained from \eqn{R} resamples. The \emph{Monte Carlo standard error} is given by 
\deqn{
\sigma(x) = \sqrt{ F(x) [1 - F(x)] / R }
}
which is the standard deviation of \eqn{G(x)}. 

A symmetric confidence interval for \eqn{F(x)} can be approximated as
\deqn{
G(x) +/- C \sigma(x)
}
where \eqn{C} is some quantile of the standard normal distribution. Note that the critical value \eqn{C} corresponds to the confidence level (\code{conf.level}) of the approximation.

Let \eqn{\alpha} denote the significance level (\code{sig.level}) for a one-sided test (\eqn{\alpha} is one-half the significance level for two-sided tests). Define \eqn{a} to be the value of the test statistic such that \eqn{F(a) = \alpha}. 

The parameter \eqn{\delta} (\code{delta}) quantifies the accuracy of the approximation, such that 
\deqn{
|G(a) - \alpha| < \alpha \delta
}
with a given confidence, which is controlled by the \code{conf.level} argument.
}
\value{
\item{mcse }{Monte Carlo standard error.}
\item{R }{Number of resamples.}
\item{delta }{Accuracy of approximation.}
\item{conf.level }{Confidence level.}
\item{sig.level }{Significance level.}
\item{alternative }{Alternative hypothesis.}
}
\references{
Helwig, N. E. (2019). Statistical nonparametric mapping: Multivariate permutation tests for location, correlation, and regression problems in neuroimaging. WIREs Computational Statistics, 11(2), e1457. doi: 10.1002/wics.1457
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
This function is only relevant for non-exact tests. For exact tests, \eqn{F(x) = G(x)} so the Monte Carlo standard error is zero.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{np.cor.test}}, \code{\link{np.loc.test}}, \code{\link{np.reg.test}}
}
\examples{

###***###   EXAMPLE 1   ###***###

# get the Monte Carlo standard error and the 
# accuracy (i.e., delta) for given R = 10000
# using the default two-sided alternative hypothesis,
# the default confidence level (conf.level = 0.95),
# and the default significance level (sig.level = 0.05)

mcse(R = 10000)

# se = 0.0016
# delta = 0.1224


###***###   EXAMPLE 2   ###***###

# get the Monte Carlo standard error and the 
# number of resamples (i.e., R) for given delta = 0.01
# using a one-sided alternative hypothesis,
# the default confidence level (conf.level = 0.95),
# and the default significance level (sig.level = 0.05)

mcse(delta = 0.1, alternative = "one.sided")

# se = 0.0026
# R = 7299

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
