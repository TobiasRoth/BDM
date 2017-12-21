\name{swissGrid}
\alias{swissGrid}

\title{
  Convert UTM or Lon/Lat coordinates to Swiss grid
}

\description{
  Convert UTM or Lon/Lat coordinates to Swiss grid.
}

\usage{
swissGrid(x, y = NA, UTM = FALSE, zone = NA, ...)
}

\arguments{
  \item{x}{Either a A \code{\link[base]{data.frame}} with columns \code{x} and \code{y} or a vector with the x coordinates.}
  \item{y}{If x is a vector, then y is a vector with the y coordinates.}
  \item{UTM}{Set \code{UTM=TRUE} if coordinates are provided as UTM coordinates, if Lon/Lat coordinates are provided then \code{UTM=FALSE}, which is the default.}
  \item{zone}{Give the \sQuote{zone} for UTM coordinates.}
  \item{\dots}{further arguments mainly \code{km=FALSE} if UTM coordinates are given in meter.}
}

\value{
A data frame with the x and y columns that contain the results of the conversion to the swiss grid.
}

\author{Tobias Roth \email{roth@hintermannweber.ch}}

\seealso{
\code{\link[PBSmapping]{convUL}}
}
\examples{
# UTM -> Swiss Grid
swissGrid(393920, 5271690, projection="UTM", zone=32, km=FALSE)
# Lon/Lat -> Swiss Grid
swissGrid(7.58901, 47.58989)
}