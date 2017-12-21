\name{Z8}
\alias{Z8}
\title{
Calculates population trends of single species
}

\description{
Calculates population trends of single species using generalized linear mixed models (GLMM) with presence/absence data or abundance data or McNemar-test.  
}

\usage{Z8(dat, surveys, occ=TRUE, abundance=FALSE, McNemar5=FALSE, McNemar10=FALSE, delNull = FALSE)}

\arguments{
  \item{dat}{The species observation data from the monitoring scheme. A \code{\link[base]{data.frame}} with five columns representing plot id, year when survey was scheduled, species id, species name and the number of individuals. The ordering of the first four columns of the \code{\link[base]{data.frame}} is mandatory. The number of individuals are not needed if the trend should be calculated only form presence/absence data. All further columns are dumped before calculation. Each row in \code{dat} represents a valid record of a species.}
  \item{surveys}{The surveys that should be analyzed. A \code{\link[base]{data.frame}} with three columns representing plot id, year when survey was scheduled and year when survey was performed. The ordering of the columns is mandatory. Each line in \code{surveys} represents a valid survey.}
  \item{occ}{Logical indicating if trend should be estimated from presence-absence data.}
  \item{abundance}{Logical indicating if trend should be estimated from abundance data.}
  \item{McNemar5}{Logical indicating if McNemar test should be used to test for significant trend within five years between the surveys on each plot.}
  \item{McNemar10}{Logical indicating if McNemar test should be used to test for significant trend within ten years between the surveys on each plot.}
  \item{delNull}{Logical indicating if surveys with no records should be removed to calculate temporal change in abundance.}
}

\details{
Using the package \code{lme4}, the function fits a generalized linear mixed model with random site and if abundance data are used also a random year effect to the data. For the presence/absence data the \sQuote{Binomial}-family and for the abundance data the \sQuote{Poisson}-familiy is used. Make sure that each line in \code{dat} contains a valid observation of a species. If the column \code{speciesID} of \code{dat}  contains \code{NAs} or is empty the respective line is dropped.
}

\value{
A list with the data.frame \code{res} and the scalar \code{N}. \code{N} gives the number of sites with two surveys in the 10 year period. In \code{res} each line gives the estimates of population trend for each recorded species. The data frame contains the following columns:
  \item{\code{speciesID}}{ID of species.}
  \item{\code{speciesName}}{Name of the species.} 
  \item{\code{Occ_1}}{Percentage of occupied sites during the first survey in the last ten years.} 
  \item{\code{Occ_2}}{Percentage of occupied sites during the second survey in the last ten years.} 
  \item{\code{N_const}}{The number of sites that were occupied during both surveys.} 
  \item{\code{N_new}}{The number of sites only occupied during the second survey.} 
  \item{\code{N_lost}}{The number of sites only occupied during the first survey.} 
  \item{\code{B-GLMM_L}}{Linear trend estimated from presence/absence data.} 
  \item{\code{B-GLMM_L_p}}{P-value of linear trend estimated from presence/absence data.} 
  \item{\code{P-GLMM_L}}{Linear trend estimated from abundance data.} 
  \item{\code{P-GLMM_L_p}}{P-value of linear trend estimated from abundance data.} 
  \item{\code{N5}}{Sample size of plots that were surveyed two times during the last ten years.} 
  \item{\code{Occ5}}{The number of occupied sites during the first survey in the last ten years.} 
  \item{\code{Increase5}}{The number of sites that were not occupied during the first survey but occupied during the second survey in the last ten years.} 
  \item{\code{Decrease5}}{The number of sites that were occupied during the first survey but not occupied during the second survey in the last ten years.} 
  \item{\code{McNemar_p5}}{P-value of trend between first and second survey in the last ten year estimated using McNemar-test.} 
 \item{\code{N10}}{Sample size of plots that were surveyed during the first and third survey during the last 15 years.} 
  \item{\code{Occ10}}{The number of occupied sites during the first survey in the last 15 years.} 
  \item{\code{Increase10}}{The number of sites that were not occupied during the first survey but occupied during the third survey in the last 15 years.} 
  \item{\code{Decrease10}}{The number of sites that were occupied during the first survey but not occupied during the third survey in the last 15 years.} 
  \item{\code{McNemar_p10}}{P-value of trend between first and third survey in the last 15 year estimated using McNemar-test.} 
}

\author{Tobias Roth \email{roth@hintermannweber.ch}}

\examples{
### Estimate population trends of Sylvia-warblers in Argovia
data(LANAG)
dat <- LANAG$bi[!is.na(LANAG$bi$speciesID) & (LANAG$bi$speciesID == 4570 | LANAG$bi$speciesID == 4600),]
surveys <- LANAG$sur[,c("coordID", "pyear", "yearBi")]
Z8(dat, surveys)

### Compare abundance and distribution trends of Coenonympha pamphilus and Lasiommata megera
data(LANAG)
dat <- LANAG$bu[!is.na(LANAG$bu$speciesID) & (LANAG$bu$speciesID == 31211 | LANAG$bu$speciesID == 31245),]
surveys <- LANAG$surBu[,c("coordID", "pyear", "yearBu")]
Z8(dat, surveys, occ=TRUE, abundance=TRUE, McNemar5=TRUE)
}