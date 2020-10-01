\name{np.boot}
\alias{np.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Bootstrap Resampling
}
\description{
Nonparametric bootstrap resampling for univariate and multivariate statistics. Computes bootstrap estimates of the standard error, bias, and covariance. Also computes five different types of bootstrap confidence intervals: normal approximation interval, basic (reverse percentile) interval, percentile interval, studentized (bootstrap-\emph{t}) interval, and bias-corrected and accelerated (BCa) interval.
}
\usage{
np.boot(x, statistic, ..., R = 9999, level = c(0.9, 0.95, 0.99),
        method = c("norm", "basic", "perc", "stud", "bca")[-4], 
        sdfun = NULL, sdrep = 99, jackknife = NULL, 
        parallel = FALSE, cl = NULL, boot.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of data (for univariate data) or vector of row indices (for multivariate data). See examples for bootstrapping multivariate data.
}
  \item{statistic}{
function that takes in \code{x} (and possibly additional arguments passed using \code{...}) and returns a vector containing the statistic(s). See examples.
}
  \item{...}{
additional named arguments for the \code{statistic} function.  
}
  \item{R}{
number of bootstrap replicates  
}
  \item{level}{
desired confidence level(s) for the computed intervals. Default computes 90\%, 95\%, and 99\% confidence intervals.  
}
  \item{method}{
method(s) for computing confidence intervals. Partial matching is allowed. Any subset of allowable methods is permitted (default computes all intervals except studentized). Set \code{method = NULL} to produce no confidence intervals.
}
  \item{sdfun}{
function for computing the standard deviation of \code{statistic}. Should produce a vector the same length as the output of \code{statistic}. Only applicable if \code{"stud" \%in\% method}. If \code{NULL}, an inner bootstrap is used to estimate the standard deviation.
}
  \item{sdrep}{
number of bootstrap replicates for the inner bootstrap used to estimate the standard deviation of \code{statistic}. Only applicable if \code{"stud" \%in\% method} and \code{sdfun = NULL}. Larger values produce more accurate estimates (see Note).
}
  \item{jackknife}{
function that takes in \code{x} (and possibly additional arguments passed using \code{...}) and returns a vector containing the jackknife statistic(s). Should produce a vector the same length as the output of \code{statistic}. Only applicable if \code{"bca" \%in\% method}. If \code{NULL}, the jackknife function is defined as the statistic function (default). See the last example for a case when \code{statistic} and \code{jackknife} are different.
}
  \item{parallel}{
Logical indicating if the \code{\link{parallel}} package should be used for parallel computing (of the bootstrap distribution). Defaults to \code{FALSE}, which implements sequential computing.
}
  \item{cl}{
Cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}. 
}
  \item{boot.dist}{
Logical indicating if the bootstrap distribution should be returned (see Note).  
}
}
\details{
The first three intervals (normal, basic, and percentile) are only first-order accurate, whereas the last two intervals (studentized and BCa) are both second-order accurate. Thus, the results from the studentized and BCa intervals tend to provide more accurate coverage rates.

Unless the standard deviation function for the studentized interval is input via the \code{sdfun} argument, the studentized interval can be quite computationally costly. This is because an inner bootstrap is needed to estimate the standard deviation of the statistic for each (outer) bootstrap replicate---and you may want to increase the default number of inner bootstrap replicates (see Note).

The efficiency of the BCa interval will depend on the sample size \eqn{n} and the computational complexity of the (jackknife) statistic estimate. Assuming that \eqn{n} is not too large and the jackknife statistic is not too difficult to compute, the BCa interval can be computed reasonably quickly---especially in comparison the studentized interval with an inner bootstrap.

Computational details of the various confidence intervals are described in Efron and Tibshirani (1994) and in Davison and Hinkley (1997). For a useful and concise discussion of the various intervals, see Carpenter and Bithell (2000).
}
\value{
\item{t0 }{Observed statistic, computed using \code{statistic(x, ...)}}
\item{se }{Bootstrap estimate of the standard error.}
\item{bias }{Bootstrap estimate of the bias.}
\item{cov }{Bootstrap estimate of the covariance (for multivariate statistics).}
\item{normal }{Normal approximation confidence interval(s).}
\item{basic }{Basic (reverse percentile) confidence interval(s).}
\item{percent }{Percentile confidence interval(s).}
\item{student }{Studentized (bootstrap-\emph{t}) confidence interval(s).}
\item{bca }{Bias-corrected and accelerated (BCa) confidence interval(s).}
\item{z0 }{Bias-correction factor(s). Only provided if \code{bca \%in\% method}.}
\item{acc }{Acceleration factor(s). Only provided if \code{bca \%in\% method}.}
\item{boot.dist }{Bootstrap distribution of statistic(s). Only provided if \code{boot.dist = TRUE}.}
\item{R }{Number of bootstrap replicates (same as input).}
\item{level }{Confidence level (same as input).}
\item{sdfun }{Standard deviation function for statistic (same as input).}
\item{sdrep }{Number of inner bootstrap replicates (same as input).}
\item{jackknife }{Jackknife function (same as input).}
}
\references{
Carpenter, J., & Bithell, J. (2000). Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians. \emph{Statistics in Medicine, 19}(9), 1141-1164. doi: 10.1002/(SICI)1097-0258(20000515)19:9\%3C1141::AID-SIM479\%3E3.0.CO;2-F 

Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and Their Application}. Cambridge University Press. doi: 10.1017/CBO9780511802843

Efron, B., & Tibshirani, R. J. (1994). \emph{An Introduction to the Boostrap}. Chapman & Hall/CRC. doi: 10.1201/9780429246593
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If \code{boot.dist = TRUE}, the output \code{boot.dist} will be a matrix of dimension \code{R} by \code{length(statistic(x, ...))} if the statistic is multivariate. Otherwise the bootstrap distribution will be a vector of length \code{R}.

For the \code{"stud"} method, the default of \code{sdrep = 99} may produce a crude estimate of the standard deviation of the statistic(s). For more accurate estimates, the value of \code{sdrep} may need to be set substantially larger, e.g., \code{sdrep = 999}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
######***######   UNIVARIATE DATA   ######***######

### Example 1: univariate statistic (median)

# generate 100 standard normal observations
set.seed(1)
n <- 100
x <- rnorm(n)

# nonparametric bootstrap
npbs <- np.boot(x = x, statistic = median)
npbs


### Example 2: multivariate statistic (quartiles)

# generate 100 standard normal observations
set.seed(1)
n <- 100
x <- rnorm(n)

# nonparametric bootstrap
npbs <- np.boot(x = x, statistic = quantile, 
                probs = c(0.25, 0.5, 0.75))
npbs



######***######   MULTIVARIATE DATA   ######***######

### Example 1: univariate statistic (correlation)

# correlation matrix square root (with rho = 0.5)
rho <- 0.5
val <- c(sqrt(1 + rho), sqrt(1 - rho))
corsqrt <- matrix(c(val[1], -val[2], val), 2, 2) / sqrt(2)

# generate 100 bivariate observations (with rho = 0.5)
n <- 100
set.seed(1)
data <- cbind(rnorm(n), rnorm(n)) \%*\% corsqrt

# define statistic function
statfun <- function(x, data) cor(data[x,1], data[x,2])

# nonparametric bootstrap
npbs <- np.boot(x = 1:n, statistic = statfun, data = data)
npbs


### Example 2: multivariate statistic (variances and covariance)

# correlation matrix square root (with rho = 0.5)
rho <- 0.5
val <- c(sqrt(1 + rho), sqrt(1 - rho))
corsqrt <- matrix(c(val[1], -val[2], val), 2, 2) / sqrt(2)

# generate 100 bivariate observations (with rho = 0.5)
n <- 100
set.seed(1)
data <- cbind(rnorm(n), rnorm(n)) \%*\% corsqrt

# define statistic function
statfun <- function(x, data) {
  cmat <- cov(data[x,])
  ltri <- lower.tri(cmat, diag = TRUE)
  cvec <- cmat[ltri]
  names(cvec) <- c("var(x1)", "cov(x1,x2)", "var(x2)")
  cvec
}

# nonparametric bootstrap
npbs <- np.boot(x = 1:n, statistic = statfun, data = data)
npbs



\dontrun{

######***######   REGRESSION   ######***######

### Example 1: bootstrap cases

# generate 100 observations
n <- 100
set.seed(1)
x <- seq(0, 1, length.out = n)
y <- 1 + 2 * x + rnorm(n)
data <- data.frame(x = x, y = y)

# define statistic function
statfun <- function(x, data) {
  xmat <- cbind(1, data$x[x])
  xinv <- solve(crossprod(xmat)) \%*\% t(xmat)
  coef <- as.numeric(xinv \%*\% data$y[x])
  names(coef) <- c("(Intercept)", "x")
  coef
}

# nonparametric bootstrap
npbs <- np.boot(x = 1:n, statistic = statfun, data = data)
npbs


### Example 2: bootstrap residuals

# generate 100 observations
n <- 100
set.seed(1)
x <- seq(0, 1, length.out = n)
y <- 1 + 2 * x + rnorm(n)

# prepare data
xmat <- cbind(1, x)
xinv <- solve(crossprod(xmat)) \%*\% t(xmat)
fit <- xmat \%*\% xinv \%*\% y
data <- list(fit = fit, resid = y - fit, xinv = xinv, x = x)

# define statistic function
statfun <- function(x, data) {
  ynew <- data$fit + data$resid[x]
  coef <- as.numeric(data$xinv \%*\% ynew)
  names(coef) <- c("(Intercept)", "x")
  coef
}

# define jackknife function
jackfun <- function(x, data){
  ynew <- data$fit[x] + data$resid[x]
  xmat <- cbind(1, data$x[x])
  xinv <- solve(crossprod(xmat)) \%*\% t(xmat)
  coef <- as.numeric(xinv \%*\% ynew)
  names(coef) <- c("(Intercept)", "x")
  coef
}

# nonparametric bootstrap
npbs <- np.boot(x = 1:n, statistic = statfun, data = data, 
                jackknife = jackfun)
npbs
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{univar}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
