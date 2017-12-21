\name{swissSHP}

\alias{swissSHP}
\alias{argovianSHP}

\docType{data}

\title{Different shapes to make distribution plots.}

\description{
Different shapes that are readily available to make distribution plots.
}
\usage{
  data(swissSHP)
  data(argovianSHP)
}

\format{
  \describe{
    \item{swissSHP}{A data frame with coordinates (Lon/Lat-coordinates) of the border of Switzerland.}
    \item{argovianSHP}{The data is an imported AGIS shapefile. It is saved as a \code{SpatialPolygonsDataFrame} of the package \code{sp}.}
  }
}

\examples{
### Plot ...
data(swissSHP)
plot(swissGrid(swissSHP), ty="l", ylab="", xlab="", las=1)

### Plot species richness of birds in Argovia
data(argovianSHP)
t <- par(mar=c(1,1,1,1))
plot(argovianSHP, col="grey93", bo="grey70", xlim=c(621000, 675000), ylim=c(223000, 273000))
data(birdsLANAG)
dat <- birdsLANAG[birdsLANAG$pyear>=2006, ]
temp <- table(dat$coordID)
points(coordID2coord(as.integer(names(temp))), pch=16, col="red", cex=temp/20)
par(t)
}

\keyword{datasets}
