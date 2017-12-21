\name{coordID2coord}
\alias{coordID2coord}

\title{Transform plot Id to coordinates of Swiss grid}
\description{
Transform plot Id to coordinates of Swiss grid.
}
\usage{
coordID2coord(coordID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordID}{Vector with format \code{xxxyyy}.}
}

\value{
A data frame with the x and y columns that contain the results of the conversion to the swiss grid.
}

\author{Tobias Roth \email{roth@hintermannweber.ch}}

\examples{
### Plot species richness of birds in Argovia
data(argovianSHP)
t <- par(mar=c(1,1,1,1))
plot(argovianSHP, col="grey93", bo="grey70", xlim=c(621000, 675000), ylim=c(223000, 273000))
data(birdsLANAG)
dat <- birdsLANAG[birdsLANAG$pyear>=2006, ]
temp <- table(dat$coordID)
points(coordID2coord(as.integer(names(temp))), pch=16, col="red", cex=temp/20)
}

