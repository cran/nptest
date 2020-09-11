\name{np.reg.test}
\alias{np.reg.test}

\title{
Nonparametric Tests of Regression Coefficients
}
\description{
Assuming a linear model of the form 
\deqn{Y = \alpha + X \beta + \epsilon}
or
\deqn{Y = \alpha + X \beta + Z \gamma + \epsilon}
this function implements permutation tests of \eqn{H_0: \beta = \beta_0} where \eqn{\beta_0} is the user-specified null vector. 
}

\usage{
np.reg.test(x, y, z = NULL, method = NULL,
            beta = NULL, homosced = FALSE,
            R = 9999, parallel = FALSE, cl = NULL,
            perm.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of predictor variables (\code{n} by \code{p}). 
}
  \item{y}{
Response vector or matrix (\code{n} by \code{m}).
}
  \item{z}{
Optional matrix of nuisance variables (\code{n} by \code{q}). 
}
  \item{method}{
Permutation method. See Details.
}
  \item{beta}{
Null hypothesis value for \eqn{\beta} (\code{p} by \code{m}). Defaults to matrix of zeros.
}
  \item{homosced}{
Are the \eqn{\epsilon} terms homoscedastic? If \code{FALSE} (default), a robust Wald test statistic is used. Otherwise the classic \eqn{F} test statistic is used.
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
With no nuisance variables in the model (i.e., \code{z = NULL}), there are three possible options for the \code{method} argument:
\tabular{lll}{
Method \tab Model \cr
perm \tab \eqn{P Y = \alpha + X \beta + \epsilon} \cr
flip \tab \eqn{S Y = \alpha + X \beta + \epsilon} \cr
both \tab \eqn{P S Y = \alpha + X \beta + \epsilon}
}
where \eqn{P} is a permutation matrix and \eqn{S} is a sign-flipping matrix.


With nuisance variables in the model, there are eight possible options for the \code{method} argument:
\tabular{lll}{
Method \tab Name \tab Model \cr
HJ \tab Huh-Jhun \tab \eqn{P Q' R_z Y = \alpha + Q' R_z X \beta + \epsilon} \cr
KC \tab Kennedy-Cade \tab \eqn{P R_z Y = \alpha + R_z X \beta + \epsilon} \cr
SW \tab Still-White \tab \eqn{P R_z Y = \alpha + X \beta + \epsilon} \cr
TB \tab ter Braak \tab \eqn{(P R_m + H_m) Y = \alpha + X \beta + Z \gamma + \epsilon} \cr
FL \tab Freedman-Lane \tab \eqn{(P R_z + H_z) Y = \alpha + X \beta + Z \gamma + \epsilon} \cr
MA \tab Manly \tab \eqn{P Y = \alpha + X \beta + Z \gamma + \epsilon} \cr
OS \tab O'Gorman-Smith \tab \eqn{Y = \alpha + P R_z X \beta + Z \gamma + \epsilon} \cr
DS \tab Draper-Stoneman \tab \eqn{Y = \alpha + P X \beta + Z \gamma + \epsilon} 
}
where \eqn{P} is permutation matrix and \eqn{Q} is defined as \eqn{R_z = Q Q'} with \eqn{Q'Q = I}. 

Note that \eqn{H_z} is the hat matrix for the nuisance variable design matrix, and \eqn{R_z = I - H_z} is the corresponding residual forming matrix. Similarly, \eqn{H_m} and \eqn{R_m} are the hat and residual forming matrices for the full model including the predictor and nuisance variables.


}

\value{
  \item{statistic }{
Test statistic value.
}
  \item{p.value }{
p-value for testing \eqn{H_0: \beta = \beta_0}.
}
  \item{perm.dist }{
Permutation distribution of \code{statistic}.
}
  \item{method }{
Permutation method.
}
  \item{null.value }{
Null hypothesis value for \eqn{\beta}.
}
  \item{homosced }{
Homoscedastic errors?
}
  \item{R }{
Number of resamples.
}
  \item{exact }{
Exact permutation test? See Note.
}
  \item{coefficients }{
Least squares estimates of \eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma} (if applicable).
}
  \item{univariate }{
Univariate test statistic value for \eqn{j}-th variable (for multivariate inputs).
}
  \item{adj.p.value }{
Adjusted p-value for testing significance of \eqn{j}-th variable (for multivariate inputs).
}
}
\references{
DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for correlation and regression coefficients. Journal of the American Statistical Association, 112(519), 1211-1220. doi: 10.1080/01621459.2016.1202117

Draper, N. R., & Stoneman, D. M. (1966). Testing for the inclusion of variables in linear regression by a randomisation technique. Technometrics, 8(4), 695-699. doi: 10.2307/1266641

Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. Journal of Business and Economic Statistics, 1(4), 292-298. doi: 10.2307/1391660

Helwig, N. E. (2019a). Statistical nonparametric mapping: Multivariate permutation tests for location, correlation, and regression problems in neuroimaging. WIREs Computational Statistics, 11(2), e1457. doi: 10.1002/wics.1457

Helwig, N. E. (2019b). Robust nonparametric tests of general linear model coefficients: A comparison of permutation methods and test statistics. NeuroImage, 201, 116030. doi: 10.1016/j.neuroimage.2019.116030

Huh, M.-H., & Jhun, M. (2001). Random permutation testing in multiple linear regression. Communications in Statistics - Theory and Methods, 30(10), 2023-2032. doi: 10.1081/STA-100106060

Kennedy, P. E., & Cade, B. S. (1996). Randomization tests for multiple regression. Communications in Statistics - Simulation and Computation, 25(4), 923-936. doi: 10.1080/03610919608813350

Manly, B. (1986). Randomization and regression methods for testing for associations with geographical, environmental and biological distances between populations. Researches on Population Ecology, 28(2), 201-218. doi: 10.1007/BF02515450

Nichols, T. E., Ridgway, G. R., Webster, M. G., & Smith, S. M. (2008). GLM permutation: nonparametric inference for arbitrary general linear models. NeuroImage, 41(S1), S72.

O'Gorman, T. W. (2005). The performance of randomization tests that use permutations of independent variables. Communications in Statistics - Simulation and Computation, 34(4), 895-908. doi: 10.1080/03610910500308230

Still, A. W., & White, A. P. (1981). The approximate randomization test as an alternative to the F test in analysis of variance. British Journal of Mathematical and Statistical Psychology, 34(2), 243-252. doi: 10.1111/j.2044-8317.1981.tb00634.x

ter Braak, C. J. F. (1992). Permutation versus bootstrap significance tests in multiple regression and ANOVA. In K. H. J\"{o}ckel, G. Rothe, & W. Sendler (Eds.), Bootstrapping and related techniques. lecture notes in economics and mathematical systems, vol 376 (pp. 79-86). Springer.

White, H. (1980). A heteroscedasticity-consistent covariance matrix and a direct test for heteroscedasticity. Econometrica, 48(4), 817-838. doi: 10.2307/1912934

Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. NeuroImage, 92, 381-397. doi: 10.1016/j.neuroimage.2014.01.060
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\note{
If \code{method = "flip"}, the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{2^n} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{2^n} possible values of the test statistic.

If \code{method = "both"}, the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{factorial(n) * (2^n)} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{factorial(n) * (2^n)} possible values of the test statistic.

If \code{method = "HJ"}, the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{factorial(n-q-1)} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{factorial(n-q-1)} possible values of the test statistic.

Otherwise the permutation test will be exact when the requested number of resamples \code{R} is greater than \code{factorial(n)} minus one. In this case, the permutation distribution \code{perm.dist} contains all \code{factorial(n)} possible values of the test statistic. 
}

\section{Multivariate Tests }{
If the input \code{y} is a matrix with \eqn{m > 1} columns, the multivariate test statistic is defined as \code{statistic = max(univariate)} given that the \code{univariate} test statistics are non-negative. 

The global null hypothesis (across all \eqn{m} variables) is tested by comparing the observed \code{statistic} to the permutation distribution \code{perm.dist}. This produces the \code{p.value} for testing the global null hypothesis.

The local null hypothesis (separately for each variable) is tested by comparing the \code{univariate} test statistic to \code{perm.dist}. This produces the adjusted p-values (\code{adj.p.values}), which control the familywise Type I error rate across the \eqn{m} tests.
}

\seealso{
\code{\link{plot.np.reg.test}} S3 plotting method for visualizing the results
}

\examples{

######******######   UNIVARIATE   ######******######

###***###   TEST ALL COEFFICIENTS   ###***###

# generate data
set.seed(1)
n <- 10
x <- cbind(rnorm(n), rnorm(n))
y <- rnorm(n)

# Wald test (method = "perm")
set.seed(0)
np.reg.test(x, y)

# F test (method = "perm")
set.seed(0)
np.reg.test(x, y, homosced = TRUE)


###***###   TEST SUBSET OF COEFFICIENTS   ###***###

# generate data
set.seed(1)
n <- 10
x <- rnorm(n)
z <- rnorm(n)
y <- 3 + 2 * z + rnorm(n)

# Wald test (method = "HJ")
set.seed(0)
np.reg.test(x, y, z)

# F test (method = "HJ")
set.seed(0)
np.reg.test(x, y, z, homosced = TRUE)


\dontrun{

######******######   MULTIVARIATE   ######******######

###***###   TEST ALL COEFFICIENTS   ###***###

# generate data
set.seed(1)
n <- 10
x <- cbind(rnorm(n), rnorm(n))
y <- matrix(rnorm(n * 3), nrow = n, ncol = 3)

# multivariate Wald test (method = "perm")
set.seed(0)
np.reg.test(x, y)

# multivariate F test (method = "perm")
set.seed(0)
np.reg.test(x, y, homosced = TRUE)


###***###   TEST SUBSET OF COEFFICIENTS   ###***###

# generate data
set.seed(1)
n <- 10
x <- rnorm(n)
z <- rnorm(n)
y <- cbind(1 + 3 * z + rnorm(n),
           2 + 2 * z + rnorm(n),
           3 + 1 * z + rnorm(n))
           
# multivariate Wald test (method = "HJ")
set.seed(0)
np.reg.test(x, y, z)

# multivariate F test (method = "HJ")
set.seed(0)
np.reg.test(x, y, z, homosced = TRUE)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ regression }
\keyword{ robust }
\keyword{ univar }
