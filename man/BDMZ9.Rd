\name{BDMZ9}
\alias{BDMZ9}
\docType{data}
\title{
Z9 data from Swiss monitoring scheme (2001-2010)
}
\description{
Biodiversity data from Swiss monitoring scheme (2001-2010) of indicator Z9 (i.e. species richness in habitats). The list includes four \code{data.frame}: one \code{data.frame} for the survey data and one for each of the three species groups (vascular plants, mooses and molluscs).
}
\usage{data(BDMZ9)}
\format{
A list with the following 4 data frames.
  \describe{
    \item{\code{sur}}{A data frame with the surveys of birds, plants and butterflies.}
    \item{\code{pla}}{A data frame with the plant observations.}
    \item{\code{moo}}{A data frame with the moos observations.}
    \item{\code{mol}}{A data frame with the mollusc observations.}
  }
}

\source{
\url{http://www.biodiversitymonitoring.ch}
}
\references{
Hintermann, U.; Weber, D.; Zangger, A., 2000: Biodiversity monitoring in Switzerland. - \emph{Schriftenreihe Landschafspflege Naturschutz} \bold{62}: 47-58.

Weber, D.; Hintermann, U.; Zangger, A., 2004: Scale and trends in species richness: considerations for monitoring biological diversity for political purposes. - \emph{Global Ecology and Biogeography} \bold{13}: 97-104.
}
\examples{
data(BDMZ9)

### Current species richness of Z9-Plants
dat <- BDMZ9$pla
surveys <- BDMZ9$sur[!is.na(BDMZ9$sur$yearPla), 1:2]
land_use <- BDMZ9$sur[!is.na(BDMZ9$sur$yearPla), "land_use"]
MAR(dat, surveys, land_use, start=2006, end=2010)

### Change in species richness of molluscs between 2001/2005 and 2006/2010
dat <- BDMZ9$mol
surveys <- BDMZ9$sur[!is.na(BDMZ9$sur$yearMol), 1:2]
land_use <- BDMZ9$sur[!is.na(BDMZ9$sur$yearMol), "land_use"]
MARdiff(dat, surveys, land_use)
}

\keyword{datasets}
