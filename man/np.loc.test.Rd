\name{np.loc.test}
\alias{np.loc.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Tests of Location Parameters
}
\description{
Performs one and two sample nonparametric (randomization) tests of location parameters, i.e., means and medians. Implements univariate and multivariate tests using eight different test statistics: Student's one-sample t-test, Johnson's modified t-test, Wilcoxon's Signed Rank test, Fisher's Sign test, Student's two-sample t-test, Welch's t-test, Wilcoxon's Rank Sum test (i.e., Mann-Whitney's \eqn{U} test), and a studentized Wilcoxon test for unequal variances.
}
\usage{
np.loc.test(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, var.equal = FALSE, 
            median.test = FALSE, symmetric = TRUE,
            R = 9999, parallel = FALSE, cl = NULL,
            perm.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector (or matrix) of data values.
}
  \item{y}{
Optional numeric vector (or matrix) of data values.
}
  \item{alternative}{
Alternative hypothesis. Must be either "two.sided" (\eqn{H_1: \mu \neq \mu_0}), "less" (\eqn{H_1: \mu < \mu_0}), or "greater" (\eqn{H_1: \mu > \mu_0}).
}
  \item{mu}{
Null hypothesis value \eqn{\mu_0}. Defaults to zero.
}
  \item{paired}{
Logical indicating whether you want a paired location test.  
}
  \item{var.equal}{
Logical indicating whether to treat the two variances as being equal.
}

  \item{median.test}{
Logical indicating whether the location test is for the median. Default is \code{FALSE}, i.e., \eqn{\mu} is the mean.
}
  \item{symmetric}{
Logical indicating if the distribution of \code{x} should be assumed to be symmetric around \eqn{\mu}. Only used for one (or paired) sample tests.
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
\tabular{ll}{
One sample \tab \eqn{\mu} is the mean (or median) of \eqn{X} \cr
Paired \tab \eqn{\mu} is the mean (or median) of \eqn{X - Y} \cr
Two sample \tab \eqn{\mu} is the mean difference \eqn{E(X) - E(Y)} \cr 
 \tab or the median of the differences \eqn{X - Y} \cr
}

For one (or paired) sample tests, the different test statistics can be obtained using
\tabular{lll}{
 \tab \code{median.test = F} \tab \code{median.test = T} \cr
\code{symmetric = F} \tab Johnson t test \tab Fisher sign test \cr
\code{symmetric = T} \tab Student t test \tab Wilcoxon signed rank test \cr
}

For two sample tests, the different test statistics can be obtained using
\tabular{lll}{
 \tab \code{median.test = F} \tab \code{median.test = T} \cr
\code{var.equal = F} \tab Welch t test \tab Studentized Wilcoxon test \cr
\code{var.equal = T} \tab Student t test \tab Wilcoxon rank sum test \cr
}


}
\value{
\item{statistic }{
Test statistic value.
}
  \item{p.value }{
p-value for testing \eqn{H_0: \mu = \mu_0}.
}
  \item{perm.dist }{
Permutation distribution of \code{statistic}.
}
  \item{alternative }{
Alternative hypothesis.
}
  \item{null.value }{
Null hypothesis value for \eqn{\mu}.
}
  \item{var.equal }{
Assuming equal variances? Only for two sample tests.
}
  \item{median.test }{
Testing the median?
}
  \item{symmetric }{
Assuming symmetry? Only for one sample and paired tests.
}
  \item{R }{
Number of resamples.
}
  \item{exact }{
Exact permutation test? See Note.
}
  \item{estimate }{
Estimate of parameter \eqn{\mu}.
}
  \item{univariate }{
Univariate test statistic value for \eqn{j}-th variable (for multivariate input).
}
  \item{adj.p.value }{
Adjusted p-value for testing significance of \eqn{j}-th variable (for multivariate input).
}
  \item{method }{
Method used for permutation test. See Details.
}
}
\references{
Blair, R. C., Higgins, J. J., Karniski, W., & Kromrey, J. D. (1994). A study of multivariate permutation tests which may replace Hotelling's T2 test in prescribed circumstances. Multivariate Behavioral Research, 29(2), 141-163. doi: 10.1207/s15327906mbr2902_2

Chung, E., & Romano, J. P. (2016). Asymptotically valid and exact permutation tests based on two-sample U-statistics. Journal of Statistical Planning and Inference, 168, 97-105. doi: 10.1016/j.jspi.2015.07.004

Fisher, R. A. (1925). Statistical methods for research workers. Edinburgh: Oliver and Boyd.

Helwig, N. E. (2019). Statistical nonparametric mapping: Multivariate permutation tests for location, correlation, and regression problems in neuroimaging. WIREs Computational Statistics, 11(2), e1457. doi: 10.1002/wics.1457

Janssen, A. (1997). Studentized permutation tests for non-i.i.d. hypotheses and the generalized Behrens-Fisher problem. Statistics & Probability Letters , 36 (1), 9-21. doi: 10.1016/S0167-7152(97)00043-6

Johnson, N. J. (1978). Modified t tests and confidence intervals for asymmetrical populations. Journal of the American Statistical Association, 73 (363), 536-544. doi: 10.2307/2286597

Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. Annals Of Mathematical Statistics, 18(1), 50-60. doi: 10.1214/aoms/1177730491

Pitman, E. J. G. (1937a). Significance tests which may be applied to samples from any populations. Supplement to the Journal of the Royal Statistical Society, 4(1), 119-130. doi: 10.2307/2984124

Romano, J. P. (1990). On the behavior of randomization tests without a group invariance assumption. Journal of the American Statistical Association, 85(411), 686-692. doi: 10.1080/01621459.1990.10474928

Student. (1908). The probable error of a mean. Biometrika, 6(1), 1-25. doi: 10.2307/2331554

Welch, B. L. (1938). The significance of the difference between two means when the population variances are unequal. Biometrika, 39(3/4), 350-362. doi: 10.2307/2332010

Wilcoxon, F. (1945). Individual comparisons by ranking methods. Biometrics Bulletin, 1(6), 80-83. doi: 10.2307/3001968
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\note{
For one sample (or paired) tests, the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{2^n} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{2^n} possible values of the test statistic.

For two sample tests, the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{choose(N, n)} minus one, where \code{m = length(x)}, \code{n = length(y)}, and \code{N = m + n}. In this case, the permutation distribution \code{perm.dist} contains all \code{choose(N, n)} possible values of the test statistic.
}

\section{Multivariate Tests }{
If the input \code{x} (and possibly \code{y}) is a matrix with \eqn{m > 1} columns, the multivariate test statistic is defined as 
\tabular{ll}{
\code{alternative} \tab \code{statistic} \cr
\code{two.sided} \tab \code{max(abs(univariate))} \cr
\code{less} \tab \code{min(univariate)} \cr
\code{greater} \tab \code{max(univariate)} \cr
}

The global null hypothesis (across all \eqn{m} variables) is tested by comparing the observed \code{statistic} to the permutation distribution \code{perm.dist}. This produces the \code{p.value} for testing the global null hypothesis.

The local null hypothesis (separately for each variable) is tested by comparing the \code{univariate} test statistic to \code{perm.dist}. This produces the adjusted p-values (\code{adj.p.values}), which control the familywise Type I error rate across the \eqn{m} tests.
}

\seealso{
\code{\link{plot.np.loc.test}} S3 plotting method for visualizing the results
}

\examples{

######******######   UNIVARIATE   ######******######

###***###   ONE SAMPLE   ###***###

# generate data
set.seed(1)
n <- 10
x <- rnorm(n, mean = 0.5)

# one sample t-test
set.seed(0)
np.loc.test(x)

# Johnson t-test
set.seed(0)
np.loc.test(x, symmetric = FALSE)

# Wilcoxon signed rank test
set.seed(0)
np.loc.test(x, median.test = TRUE)

# Fisher sign test
set.seed(0)
np.loc.test(x, median.test = TRUE, symmetric = FALSE)


###***###   PAIRED SAMPLE   ###***###

# generate data
set.seed(1)
n <- 10
x <- rnorm(n, mean = 0.5)
y <- rnorm(n)

# paired t-test
set.seed(0)
np.loc.test(x, y, paired = TRUE)

# paired Johnson t-test
set.seed(0)
np.loc.test(x, y, paired = TRUE, symmetric = FALSE)

# paired Wilcoxon signed rank test
set.seed(0)
np.loc.test(x, y, paired = TRUE, median.test = TRUE)

# paired Fisher sign test
set.seed(0)
np.loc.test(x, y, paired = TRUE, median.test = TRUE, symmetric = FALSE)


###***###   TWO SAMPLE   ###***###

# generate data
set.seed(1)
m <- 7
n <- 8
x <- rnorm(m, mean = 0.5)
y <- rnorm(n)

# Welch t-test
set.seed(0)
np.loc.test(x, y)

# Student t-test
set.seed(0)
np.loc.test(x, y, var.equal = TRUE)

# Studentized Wilcoxon test
set.seed(0)
np.loc.test(x, y, median.test = TRUE)

# Wilcoxon rank sum test
set.seed(0)
np.loc.test(x, y, var.equal = TRUE, median.test = TRUE)



\dontrun{

######******######   MULTIVARIATE   ######******######

###***###   ONE SAMPLE   ###***###

# generate data
set.seed(1)
n <- 10
x <- cbind(rnorm(n, mean = 0.5), 
           rnorm(n, mean = 1), 
           rnorm(n, mean = 1.5))

# multivariate one sample t-test
set.seed(0)
ptest <- np.loc.test(x)
ptest
ptest$univariate
ptest$adj.p.values


###***###   PAIRED SAMPLE   ###***###

# generate data
set.seed(1)
n <- 10
x <- cbind(rnorm(n, mean = 0.5), 
           rnorm(n, mean = 1), 
           rnorm(n, mean = 1.5))
y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)

# multivariate paired t-test
set.seed(0)
ptest <- np.loc.test(x, y, paired = TRUE)
ptest
ptest$univariate
ptest$adj.p.values


###***###   TWO SAMPLE   ###***###

# generate data
set.seed(1)
m <- 7
n <- 8
x <- cbind(rnorm(m, mean = 0.5), 
           rnorm(m, mean = 1), 
           rnorm(m, mean = 1.5))
y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)

# multivariate Welch t-test
set.seed(0)
ptest <- np.loc.test(x, y)
ptest
ptest$univariate
ptest$adj.p.values

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ robust }
\keyword{ univar }
