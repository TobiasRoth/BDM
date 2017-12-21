\name{MARdiff}
\alias{MARdiff}

\title{Calculates differences in mean species richness between two periods}

\description{
  Calculates differences in mean species richness between two periods
  Calculates differences in mean species richness between two periods for a selection of plots and a given time interval. The function may also calculate the mean species richness for a selection of species (i.e. different species groups) as well as a different groups of plots (i.e. plots with different land use).
}
\usage{
MARdiff(dat, surveys, land_use = NA, species_group = NA, 
  P1 = c(2001, 2005), P2 = c(2006, 2010))
}


\arguments{
  \item{dat}{The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with three columns representing plot id, year and species id. The ordering of the first three columns of the \code{\link[base]{data.frame}} is mandatory. All further columns are dumped before calculation. Each row in \code{dat} represents the valid record of a species.}
  \item{surveys}{The surveys that should be analyzed. A \code{\link[base]{data.frame}} with two columns representing plot id and year of all the surveys that should be analyzed (including surveys with no recorded species). The ordering of the columns is mandatory. Each line in \code{surveys} represents a valid survey.}
  \item{land_use}{A vector or factor with length equal to \code{nrow(surveys)} that groups the surveys to different groups. It may be land use types or biogeographical regions.}
  \item{species_group}{A \code{\link[base]{data.frame}} with two columns representing species id and species group of all the species. The ordering of the columns is mandatory. Each row in \code{species_group} represents a species.}
  \item{P1}{The first and the last year of the first period.}
  \item{P2}{The first and the last year of the second period.}
}

\details{Make sure that each line in \code{dat} contains a valid observation of a species. If an individual is not identified on species level, then the column \code{speciesID} of \code{dat} may contain \code{NAs}.
The \code{surveys} makes sure that surveys with no recorded species are considered in the calculation of the mean species richness as well.}

\value{The function returns a list with
  \item{res}{A \code{\link[base]{data.frame}} with the results of mean species richness (\code{MAR}), lower (\code{CIL}) and upper (\code{CIU}) 95\%-confidence interval for each 5-year period and each strata. The column \code{N_species} contains the number of recorded species. If calculated, a column for land use and/or species group are also provided.}
  \item{land_use_change}{A \code{\link[base]{table}} with the number of plots per land use type during the firs and second period. Linest give land use types for the first periods and columns give the land use type for the second period.}
}

\author{Tobias Roth \email{roth@hintermannweber.ch}}

\examples{
## calculate mean species richness over time
data(birdsLANAG)
dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
surveys$pyear <- substr(surveys$coordID, 8,12)
surveys$coordID <- substr(surveys$coordID, 1,6)
MARdiff(dat, surveys)
}