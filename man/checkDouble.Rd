\name{checkDouble}
\alias{checkDouble}
\title{Check for species that are recorded two times during a survey}
\description{Check whether one species is recorded several times during one survey}
\usage{checkDouble(dat)}

\arguments{
  \item{dat}{The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with three columns representing plot id, year and species id.The ordering of the first three row of the \code{\link[base]{data.frame}} is mandatory. All further columns are dumped before calculation. Each line in the \code{dat} represents the valid record of a species.}
}
\details{Make sure that each line in \code{dat} contains a valid observation of a species. If an individual is not identified on species level, then the column \code{speciesID} of \code{dat} may contain \code{NAs} Observations that are not identified on species level will be droped before check.
}

\value{Gives only a message if there were no species recorded two times during a survey. If there where surveys with species that were recorded severeal times, the function returns a vector with the survey ids of the surveys where species have recorded several times. The survey id has the following format \code{xxxyyy_year}, where \code{xxxyyy} is the plot id.
}

\author{  Tobias Roth \email{roth@hintermannweber.ch}}

\examples{
## calculate mean species richness over time
data(birdsLANAG)
checkDouble(birdsLANAG[,c("coordID", "pyear", "speciesID")])
}

