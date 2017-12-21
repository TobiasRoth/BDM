\name{LANAG}
\alias{LANAG}
\docType{data}
\title{
Biodiversity data from Agrovian monitoring scheme (1995-2011)
}
\description{
Biodiversity data from Agrovian monitoring scheme (1995-2011). The list includes six \code{data.frame}: two \code{data.frame} for the survey data of birds/plants/molluscs and butterflies, respectively and one for each of the four species groups (birds, plants, molluscs and butterflies).
}
\usage{data(LANAG)}
\format{
A list with the following 6 data frames.
  \describe{
    \item{\code{sur}}{A data frame with the surveys of birds, plants and butterflies.}
    \item{\code{surBu}}{A data frame with the surveys of butterflies.}
    \item{\code{bi}}{A data frame with the bird observations.}
    \item{\code{pl}}{A data frame with the plasnt observations.}
    \item{\code{mo}}{A data frame with the mollusc observations.}
    \item{\code{bu}}{A data frame with the butterfly observations.}
  }
}

\source{
\url{http://www.ag.ch/alg/de/pub/natur_landschaft/erfolgskontrolle/lanag.php}
}
\references{
Roth, T.; Amrhein, V.; Peter, B.; Weber, D., 2008: A Swiss agri-environment scheme effectively enhances species richness for some taxa over time. \emph{Agriculture Ecosystems & Environment} \bold{125}: 167-172.
}
\examples{
data(LANAG)
summaryBDM(LANAG$bi)
}
\keyword{datasets}
