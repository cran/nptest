\name{plot}
\alias{plot}
\alias{plot.np.cor.test}
\alias{plot.np.loc.test}
\alias{plot.np.reg.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots Permutation Distribution for Nonparametric Tests
}
\description{
\code{plot} methods for object classes "np.cor.test", "np.loc.test", and "np.reg.test"
}
\usage{
\method{plot}{np.cor.test}(x, alpha = 0.05, col = "grey", col.rr = "red",
     col.stat = "black", lty.stat = 2, lwd.stat = 2, 
     xlab = "Test Statistic", main = "Permutation Distribution", 
     breaks = "scott", border = NA, box = TRUE, ...)
           
\method{plot}{np.loc.test}(x, alpha = 0.05, col = "grey", col.rr = "red",
     col.stat = "black", lty.stat = 2, lwd.stat = 2, 
     xlab = "Test Statistic", main = "Permutation Distribution", 
     breaks = "scott", border = NA, box = TRUE, ...)
           
\method{plot}{np.reg.test}(x, alpha = 0.05, col = "grey", col.rr = "red",
     col.stat = "black", lty.stat = 2, lwd.stat = 2, 
     xlab = "Test Statistic", main = "Permutation Distribution", 
     breaks = "scott", border = NA, box = TRUE, SQRT = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "np.cor.test" output by the \code{\link{np.cor.test}} function, "np.loc.test" output by the \code{\link{np.loc.test}} function, or "np.reg.test" output by the \code{\link{np.reg.test}} function
}
  \item{alpha}{
significance level of the nonparametric test
}
  \item{col}{
color for plotting the non-rejection region
}
  \item{col.rr}{
color for plotting the rejection region
}
  \item{col.stat}{
color for plotting the observed test statistic
}
  \item{lty.stat}{
line type for plotting the observed test statistic
}
  \item{lwd.stat}{
line width for plotting the observed test statistic
}
  \item{xlab}{
x-axis label for the plot
}
  \item{main}{
title for the plot
}
  \item{breaks}{
defines the breaks of the histogram (see \code{\link{hist}})
}
  \item{border}{
color of the border around the bars
}
  \item{box}{
should a box be drawn around the plot?
}
  \item{SQRT}{
for regression tests, should the permutation distribution (and test statistic) be plotted on the square-root scale?
}
  \item{...}{
additional arguments to be passed to \code{\link{hist}}  
}
}
\details{
Plots a histogram of the permutation distribution and the observed test statistic. The argument 'alpha' controls the rejection region of the nonparametric test, which is plotted using a separate color (default is red).
}
\references{
Helwig, N. E. (2019). Statistical nonparametric mapping: Multivariate permutation tests for location, correlation, and regression problems in neuroimaging. WIREs Computational Statistics, 11(2), e1457. doi: 10.1002/wics.1457
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{np.cor.test}} for information on nonparametric correlation tests

\code{\link{np.loc.test}} for information on nonparametric location tests

\code{\link{np.reg.test}} for information on nonparametric regression tests
}
\examples{

######******######   np.cor.test   ######******######

# generate data
rho <- 0.5
val <- c(sqrt(1 + rho), sqrt(1 - rho))
corsqrt <- matrix(c(val[1], -val[2], val), 2, 2) / sqrt(2)
set.seed(1)
n <- 50
z <- cbind(rnorm(n), rnorm(n)) \%*\% corsqrt
x <- z[,1]
y <- z[,2]

# test H0: rho = 0
set.seed(0)
test <- np.cor.test(x, y)

# plot results
plot(test)


######******######   np.loc.test   ######******######

# generate data
set.seed(1)
n <- 50
x <- rnorm(n, mean = 0.5)

# one sample t-test
set.seed(0)
test <- np.loc.test(x)

# plot results
plot(test)


######******######   np.reg.test   ######******######

# generate data
set.seed(1)
n <- 50
x <- cbind(rnorm(n), rnorm(n))
beta <- c(0.25, 0.5)
y <- x \%*\% beta + rnorm(n)

# Wald test (method = "perm")
set.seed(0)
test <- np.reg.test(x, y)

# plot results
plot(test)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ htest }
\keyword{ nonparametric }
\keyword{ aplot }
\keyword{ dplot }
