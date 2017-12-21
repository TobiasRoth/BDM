\name{MAR}
\alias{MAR}
\title{Calculates mean species richness}

\description{Calculates mean species richness for for a selection of plots and a given time interval. The function may also calculate the mean species richness for a selection of species (i.e. different species groups) as well as a selection of plots (i.e. plots with different land use). }

\usage{MAR(dat, surveys, land_use=NA, species_group=NA, 
	start=2001, end=2010, period=5)}
 
\arguments{
  \item{dat}{The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with three columns representing plot id, year and species id. The ordering of the first three columns of the \code{\link[base]{data.frame}} is mandatory. All further columns are dumped before calculation. Each row in \code{dat} represents the valid record of a species.}
  \item{surveys}{The surveys that should be analyzed. A \code{\link[base]{data.frame}} with two columns representing plot id and year of all the surveys that should be analyzed (including surveys with no recorded species). The ordering of the columns is mandatory. Each line in \code{surveys} represents a valid survey.}
  \item{land_use}{A vector or factor with length equal to \code{nrow(surveys)} that groups the surveys to different groups. It may be land use types or biogeographical regions.}
  \item{species_group}{A \code{\link[base]{data.frame}} with two columns representing species id and species group of all the species. The ordering of the columns is mandatory. Each row in \code{species_group} represents a species.}
  \item{start}{The starting year from the first period the mean species richness should be calculated.}
  \item{end}{The last year of the last period the mean species richness should be calculated.}
  \item{period}{The time period for which the mean species richness should be calculated.}
}
\details{
Make sure that each line in \code{dat} contains a valid observation of a species. If an individual is not identified on species level, then the column \code{speciesID} of \code{dat} may contain \code{NAs}.

The \code{surveys} makes sure that surveys with no recorded species are considered in the calculation of the mean species richness as well.

Be aware that when mean species richness are calculated for species groups unidentified species are droped before   calculation.
}

\value{
  \item{res}{A \code{\link[base]{data.frame}} with the results of mean species richness (\code{MAR}), lower (\code{CIL}) and upper (\code{CIU}) 95\%-confidence interval for each 5-year period and each strata. The column \code{REF} contains the mean of the 25\% best plots. The column \code{N_species} contains the number of recorded species. If calculated, a column for land use and/or species group are also provided.}
}

\author{  Tobias Roth \email{roth@hintermannweber.ch}}
\examples{
## calculate mean species richness over time
data(birdsLANAG)
dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
surveys$pyear <- substr(surveys$coordID, 8,12)
surveys$coordID <- substr(surveys$coordID, 1,6)
MAR(dat, surveys)
}

