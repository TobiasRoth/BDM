\name{birdsLANAG}
\alias{birdsLANAG}
\docType{data}
\title{Birds data from Agrovian monitoring scheme (1995-2010)}
\description{Data from monitoring scheme of Agrovia.}
\usage{data(birdsLANAG)}
\format{
  A data frame with 30693 observations on the following 8 variables.
  \describe{
    \item{\code{coordID}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{pyear}}{a numeric vector}
    \item{\code{speciesID}}{a numeric vector}
    \item{\code{speciesName}}{a factor}
    \item{\code{land_use}}{a factor with levels 
          \code{keine_Hauptnutzung} 
          \code{Landwirtschaft} 
          \code{Siedlung} 
          \code{Wald}}
    \item{\code{altitude}}{a numeric vector}
    \item{\code{BGR}}{a factor with levels 
          \code{Flusstal } 
          \code{H?gelland } 
          \code{Jura }}
  }
}
\examples{
## calculate mean species richness over time
## calculate mean species richness over time
data(birdsLANAG)
dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
surveys$pyear <- substr(surveys$coordID, 8,12)
surveys$coordID <- substr(surveys$coordID, 1,6)
MAR(dat, surveys)
}
\keyword{datasets}
