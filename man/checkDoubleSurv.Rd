\name{checkDoubleSurv}
\alias{checkDoubleSurv}
\title{Check whether plots are surveyed more than once whithin any 5-year interval}
\description{Check whether plots are surveyed more than once within any 5-year interval}
\usage{checkDoubleSurv(dat)}

\arguments{
  \item{dat}{The survey data (Kopfdaten) from the monitoring scheme. A \code{\link[base]{data.frame}} with two columns representing plot id and year.The ordering of the rows of the \code{\link[base]{data.frame}} is mandatory. All further columns are dumped before calculation. Each line in the \code{dat} represents a valid survey.}
}
\details{Each plot is planned to be surveyed once within an interval of five years. checkDoubleSurv checks whether there are any plots in \code{dat} that were surveyed more than once within any interval of five years according to Planungsjahr.
}

\value{If at least one plot is surveyed more than once within any interval of five years, the output of the function is a data frame containing CoordID and Planungsjahr of the survey. If no plot was surveyed more than once within any interval of five years, the number of annual surveys is listed.
}

\author{  Nicolas Strebel \email{strebel@hintermannweber.ch}}

\examples{
## Check whether BDM Z9 plots are surveyed more than once in any 5-year interval
data(BDMZ9)
checkDoubleSurv(BDMZ9$sur[,c("coordID", "pyear")])
}

