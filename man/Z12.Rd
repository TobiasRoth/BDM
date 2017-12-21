\name{Z12}
\alias{Z12}
\title{Calculates species dissimilarity between a selection of plots}
\description{Calculates species dissimilarity from a given stratum and a given time interval. The function may be used to calculate species dissimilarity of different groups of plots (e.g. plots with different land use) or for a selection of species (i.e. different species groups).}
\usage{Z12(dat, surveys, land_use=NA, species_group=NA, start=2001, end=2010, 
  N_sim=0, signif=2, period=5, method="simpson")}

\arguments{
  \item{dat}{The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with three columns representing plot id, year and species id. The ordering of the first three columns of the \code{\link[base]{data.frame}} is mandatory. All further columns are dumped before calculation. Each row in \code{dat} represents the valid record of a species.}
  \item{surveys}{The surveys that should be analyzed. A \code{\link[base]{data.frame}} with two columns representing plot id and year of all the surveys that should be analyzed (surveys with no recorded species will not be used). The ordering of the columns is mandatory. Each line in \code{surveys} represents a valid survey.}
  \item{land_use}{A vector or factor with length equal to \code{nrow(surveys)} that groups the surveys to different groups. It may be land use types or biogeographical regions.}
  \item{species_group}{A \code{\link[base]{data.frame}} with two columns representing species id and species group of all the species. The ordering of the columns is mandatory. Each row represents a species.}
  \item{start}{The starting year from the first period the species dissimilarity should be calculated.}
  \item{end}{The last year of the last period the species dissimilarity should be calculated.}
  \item{N_sim}{The number of simulations that should for the simulation routine to account for not identified observations.}
  \item{signif}{Integer indicating the number of decimal places to be used for the index value.}
  \item{period}{The time period for which the species dissimilarity should be calculated.}
  \item{method}{ Binary Similarity index that is used for the calculation. The argument is from the function \code{\link[simba]{sim}} of the package \code{simba}.}
}

\details{
Make sure that each line in \code{dat} contains a valid observation of a species. If an individual is not identified on species level, then the column \code{speciesID} of \code{dat} may contain \code{NAs}.

The \code{surveys} makes sure that surveys with no recorded species are considered in the calculation of the mean species richness as well.

Be aware that when species dissimilarity are calculated for species groups unidentified species are droped before the  calculation of dissimilarity index. 
}

\value{
  \item{res}{A \code{\link[base]{data.frame}} with the results of mean species richness (\code{MAR}), lower (\code{CIL}) and upper (\code{CIU}) 95\%-confidence interval for each 5-year period and each strata. The column \code{N_species} contains the number of recorded species. If calculated, a column for land use and/or species group are also provided.}
}

\author{  Tobias Roth \email{roth@hintermannweber.ch}}

\references{
  Buhler, C., and T. Roth. 2011. Spread of common species results in local-scale floristic homogenization in grassland of Switzerland. \emph{Diversity and Distributions} \bold{17}:1089-1098.
}

\note{ The similarity/dissimilarity index is calculated using the function \code{\link[simba]{sim}} of the package \code{simba}.
}

\seealso{ 
  \code{\link[simba]{sim}}}

\examples{
## calculate mean species richness over time
data(birdsLANAG)
dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
surveys$pyear <- substr(surveys$coordID, 8,12)
surveys$coordID <- substr(surveys$coordID, 1,6)
Z12(dat, surveys, start=2005, end=2010)
}

