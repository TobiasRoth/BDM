\name{aggrButterflies}
\alias{aggrButterflies}

\title{
Aggregates butterfly data from Swiss BDM
}
\description{
Aggregates butterfly data from Swiss BDM to valid taxa.
}

\usage{
aggrButterflies(dat)
}

\arguments{
  \item{dat}{Vector of the species IDs to be aggregated.}
}

\value{
  \item{dat}{Vector of the species IDs aggregated to the relevant species}
}

\examples{
d <- c(31256, 31013, 31011)
d <- aggrButterflies(d)
d
}
