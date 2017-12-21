\name{BDM-package}
\alias{BDM-package}
\alias{BDM}
\docType{package}
\title{Collection of functions to analyze biodiversity data from Switzerland}
\description{	Collection of functions to analyze data from different biodiversity monitoring programs from Switzerland. 
		The package includes functions to calculate baseline data, functions for data controlling, 
		functions for data handling and functions to make figures. Furthermore it includes some real data from
		different monitoring schemes in Switzerland. Including BDM Data until 2012.}
\details{
	\tabular{ll}{
		Package: \tab BDM\cr
		Type: \tab Package\cr
		Version: \tab 0.1-7\cr
		Date: \tab 2013-02-12\cr
		License: \tab GPL (>= 2)
		}
	}

\author{	Tobias Roth \email{roth@hintermannweber.ch}
	}

\references{	\url{www.biodiversitymonitoring.ch} and
		\url{www.hintermannweber.ch}
	}
\examples{
## calculate mean species richness over time
data(birdsLANAG)
dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
surveys$pyear <- substr(surveys$coordID, 8,12)
surveys$coordID <- substr(surveys$coordID, 1,6)
MAR(dat, surveys)
	}
\keyword{ package }