\name{logit}
\alias{logit}
\alias{probit}

\title{
Logit- and probit functions
}
\description{
Function for probit or logit transformation.
}

\usage{
logit(x)
probit(x)
}

\arguments{
  \item{x}{Scalar or vector of the values that should be transformed.}
}

\value{
  \item{log(x/(1-x))}{logit transformation}
  \item{exp(x)/(1+exp(x))}{probit transformation}
}

\examples{
logit(0.3)
probit(-0.8472979)
}
