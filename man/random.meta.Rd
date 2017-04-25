\name{random.meta}
\alias{random.meta}

\title{
Exact Inference for Meta Analysis With Random Effects Model
}

\description{
Computes the point estimator for the center (theta), the p value for testing if the center is zero, and the 95\% confidence interval in a random effects model meta analysis. When the number of studies is moderate or small (<=20),  the exact inference results are based on the exact computation. When the number of studies is big (>20), the exact inference results are based on Monte-Carlo simulation.
}


\usage{
random.meta(y, v, type="DL", B=500, N=10000, Bstep=5, plot.meta=T)
}


\arguments{
  \item{y}{
A vector of the respective estimators of the study-specific effect from each study. Length should be the same as the number of studies.
}

  \item{v}{
A vector with the variance of each estimator in y. Length should be the same as the number of studies.
}

  \item{type}{The test method to be used for constructing the CI, choosing from "DL", "wang", "median" and "wilcox".  The default is "DL".
}

 \item{B}{The number of grids used to construct the 95\% CI. The default value is 500.
}

 \item{N}{The number of simulations in the Monte-Carlo simulation. The default value is 10000.
}

 \item{Bstep}{The number of steps used in searching the endpoint of the 95\% CI. The default value is 5, which the user does not need to adjust. A larger value may slow down computation.
}

 \item{plot.meta}{The logic value for generating the forest plot of the meta analysis. The default value is "TRUE".

}

}
\details{
The inference results are "exact" if K <= 20 and based on Monte-Carlo simulation if K > 20.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{theta }{The point estimator for the center}
 \item{pvalue }{The p value for testing if the center is zero}
 \item{ci95 }{The 95\% CI for the center}
}
\references{
Sifan Liu, Lu Tian, Steve Lee and Min-ge Xie (2016)  Exact inference on meta-analysis with generalized fixed-effects and random-effects models.  Tech Report.\\

Yan Wang and  Lu Tian (2017)  An efficient numerical algorithm for exact inference in meta analysis. Tech Report.
}
\author{
Lu Tian and Grace Deng
}

\examples{
##### Generate the data for a meta analysis with 8 studies #####
set.seed(100)
K=8
tau=2
v=rchisq(K, 3)
y=rnorm(K)*sqrt(v+tau)+1

##### Exact inference using the DL method #####
fit=random.meta(y, v, type="DL")
fit

##### Exact inference using the Wilcoxon method #####
fit=random.meta(y, v, type="wilcox")
fit



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Exact Inference}
\keyword{ Meta Analysis}
\keyword{ Random Effect}
