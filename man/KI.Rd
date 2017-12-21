\name{KI}
\alias{KI}
\title{Aggregates the species richness of species groups to an index}

\description{Aggregates the species richness of several species groups to a single index. The function divides the mean species richness by the value provided in 'base' and takes the geometric mean of all species groups. If the species richness equals the values of 'base' then the calculated index is 100.}

\usage{KI(dat, base, signif=0, Z12=TRUE)}

\arguments{
  \item{dat}{A list with an element for each species group that the index should be based on. Each element ist the output of the function \code{MAR()}.}
  \item{base}{Vector with the base-value for each species group.}
  \item{signif}{Integer indicating the number of decimal places to be used for the index value.}
  \item{Z12}{Logical indicating wether the index is based on Z12 values (TRUE) or on MAR values (FALSE).}
}
\details{
Note that the names of \code{base} should be the same as the names of the elements in \code{dat}, thus \code{names(dat)} should be exactly be the same as \code{base}.
}

\value{
  \item{res}{A \code{\link[base]{data.frame}} with the base value for each species group (\code{MARx}) and the aggregated species richness over all species groups (\code{KI}).}
}

\author{  Tobias Roth \email{roth@hintermannweber.ch}}
\examples{
data(LANAG)
sur <- LANAG$sur
surBU <- LANAG$surBu
bi <- LANAG$bi
pl <- LANAG$pl
mo <- LANAG$mo
bu <- LANAG$bu
resmar <- list(
  bi=MAR(bi, sur[!is.na(sur$yearBi),1:2], sur[!is.na(sur$yearBi),"land_use_bi"], start=1996, end=2011),
  pl=MAR(pl, sur[!is.na(sur$yearPf),1:2], sur[!is.na(sur$yearPf),"land_use"], start=1996, end=2011),
  mo=MAR(mo, sur[!is.na(sur$yearMo),1:2], sur[!is.na(sur$yearMo),"land_use"], start=1996, end=2011),
  bu=MAR(bu, surBU[!is.na(surBU$yearBu),1:2], surBU[!is.na(surBU$yearBu),"land_use"], start=1996, end=2011)
)
resmar$bi <- resmar$bi[resmar$bi$land_use!="keine_Hauptnutzung",]
resmar$bu <- resmar$bu[resmar$bu$land_use!="keine_Hauptnutzung",]
base <- c(11.7, 14.0, 6.6, 6.1)
names(base) <- c("bi", "pl", "mo", "bu")
KI(resmar, base)
}

