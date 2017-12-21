\name{summaryBDM}
\alias{summaryBDM}
\title{Function that gives some useful summary statistics of BDM data}
\description{Function that returns some useful summary statistics of biodiversity data, such as the number of survey, number of plots, species richness for each year, percentage of unidentified species, the number of records of each species and more.}
\usage{summaryBDM(dat)}
\arguments{
  \item{dat}{
    The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with four columns representing plot id, year, species id and species name.The ordering of the first four row of the data frame is mandatory. All further columns are dumped before calculation. Each line in the \code{dat} represents the valid record of a species.
  }
}

\value{
  Returns a list with the following objects:
    \item{n_surveys }{The number of surveys that were conducted (with at least one species recorded).}
    \item{n_plots }{The number of plots on which the \code{n_surveys} were conducted.}
    \item{SR_mean }{The yearly mean numbers of species recorded per survey.}
    \item{SR_sd }{The yearly standard deviation of the numbers of species recorded per survey.}
    \item{n_rec_species }{The number of species recorded totally.}
    \item{unident_species }{The percentages of unidentified species over all surveys (mean and standard deviation are given).}
    \item{species }{A \code{\link[base]{data.frame}} that gives for each species the number of records.}

}

\note{
  Note that surveys with no recorded species are not included into \code{dat}. Therefore, these surveys are also not accounted for in the summary statistics of \code{summaryBDM(dat)}.
}

\author{
  Tobias Roth \email{roth@hintermannweber.ch}
}

\examples{
## calculates some useful summary statistics of a BDM data-set
data(birdsLANAG)
summaryBDM(birdsLANAG[,c("coordID", "pyear", "speciesID", "speciesName")])
}
