#'Calculates mean species richness
#'
#'Calculates mean species richness for for a selection of plots and a given time
#'interval. The function may also calculate the mean species richness for a
#'selection of species (i.e. different species groups) as well as a selection of
#'plots (i.e. plots with different land use).
#'
#'Make sure that each line in \code{dat} contains a valid observation of a
#'species. If an individual is not identified on species level, then the column
#'\code{speciesID} of \code{dat} may contain \code{NAs}. The \code{surveys}
#'makes sure that surveys with no recorded species are considered in the
#'calculation of the mean species richness as well. Be aware that when mean
#'species richness are calculated for species groups unidentified species are
#'droped before   calculation.
#'
#'@param dat The species observation data from the monitoring scheme. A
#'  \code{\link[base]{data.frame}} with three columns representing plot id, year
#'  and species id. The ordering of the first three columns of the
#'  \code{\link[base]{data.frame}} is mandatory. All further columns are dumped
#'  before calculation. Each row in \code{dat} represents the valid record of a
#'  species.
#'@param surveys The surveys that should be analyzed. A
#'  \code{\link[base]{data.frame}} with two columns representing plot id and
#'  year of all the surveys that should be analyzed (including surveys with no
#'  recorded species). The ordering of the columns is mandatory. Each line in
#'  \code{surveys} represents a valid survey.
#'@param land_use A vector or factor with length equal to \code{nrow(surveys)}
#'  that groups the surveys to different groups. It may be land use types or
#'  biogeographical regions.
#'@param species_group A \code{\link[base]{data.frame}} with two columns
#'  representing species id and species group of all the species. The ordering
#'  of the columns is mandatory. Each row in \code{species_group} represents a
#'  species.
#'@param start The starting year from the first period the mean species richness
#'  should be calculated.
#'@param end The last year of the last period the mean species richness should
#'  be calculated.
#'@param period The time period for which the mean species richness should be
#'  calculated.
#'
#'@return {A \code{\link[base]{data.frame}} with the results of mean species
#'  richness (\code{MAR}), lower (\code{CIL}) and upper (\code{CIU})
#'  95\%-confidence interval for each 5-year period and each strata. The column
#'  \code{REF} contains the mean of the 25\% best plots. The column
#'  \code{N_species} contains the number of recorded species. If calculated, a
#'  column for land use and/or species group are also provided.}
#'
#' @examples
#' # Calculate mean species richness over time
#' data(birdsLANAG)
#' dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
#' surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
#' surveys$pyear <- substr(surveys$coordID, 8,12)
#' surveys$coordID <- substr(surveys$coordID, 1,6)
#' MAR(dat, surveys)
#' 
#'@export 
#' 
MAR <- function(dat, surveys, land_use=NA, species_group=NA, start=2001, end=2010, period=5) {
  
  ### Select relevant data from 'dat' and make unique ID for survey	
  if(class(dat) != "data.frame") stop("'dat' should be a data.frame")
  if(ncol(dat)<3) stop("'dat' should contain three columns")
  dat <- dat[,1:3]
  names(dat) <- c("coordID", "pyear", "speciesID")
  dat$surveyID <- paste(dat$coordID, "_", dat$pyear, sep="")
  
  ### Select relevant data from 'surveys' and make unique ID for survey
  if(class(surveys) != "data.frame") stop("'surveys' should be a data.frame")
  if(ncol(surveys)!=2) stop("'surveys' should contain the two columns 'coordID' and 'pyear'")
  names(surveys) <- c("coordID", "pyear")
  surveys$surveyID <- paste(surveys$coordID, "_", surveys$pyear, sep="")
  if(sum(table(surveys$surveyID)>1)>0) stop("There are more than one entry for a survey. Check in 'survey'")
  row.names(surveys) <- surveys$surveyID
  
  ### Prepare 'land_use'
  if(!is.na(land_use)[1]) {
    surveys$land_use <- land_use
    dat$land_use <- surveys[as.character(dat$surveyID),"land_use"]
  }
  
  ### Prepare 'species_group'
  if(!is.na(species_group)[1]) {
    if(!class(species_group)=="data.frame") stop("'species_group' should be a data.frame")
    if(ncol(species_group)!=2) stop("'species_group' should contain the two columns 'coordID' and 'group'")
    names(species_group) <- c("speciesID", "group")
    if(sum(table(species_group$speciesID)>1)>0) stop("Some species in 'species_group' are mentioned twice." )
    row.names(species_group) <- as.character(species_group$speciesID)
    dat$group <- species_group[as.character(dat$speciesID),"group"]
  }
  
  ### Prepare data.frame 'res' for the results and set counter 'i' to one
  res <- data.frame()
  i <- 1
  
  ### Calculate MAR for all records and every land use
  for(t in start:(end-(period-1))) {
    d  <- surveys[surveys$pyear >= t & surveys$pyear <= (t+(period-1)), ]
    d$SR <- 0
    temp <- table(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)), "surveyID"])
    if(nrow(d) < length(temp)) message("The data.frame 'dat' contains observations from surveys that are not contained in 'surveys'.")
    d[names(temp),"SR"] <- temp
    d <- d[!is.na(d$coordID),]
    res[i, "start"] <- t
    res[i, "end"] <- t+(period-1)
    res[i,"land_use"] <- "All"
    if(!is.na(species_group)[1]) res[i,"spec_group"] <- "All"
    res[i, "N"] <- nrow(d)
    res[i, "MAR"] <- round(mean(d$SR), 2)
    try({
      res[i, "CIL"] <- round(t.test(d$SR)$conf.int[1], 2)
      res[i, "CIU"] <- round(t.test(d$SR)$conf.int[2], 2)
    })
    res[i, "REF"] <- round(mean(d$SR[d$SR >= quantile(d$SR, 0.75)]), 2)
    speciesID <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)), "speciesID"]
    res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
    i <- i+1
    if(!is.na(land_use)[1]) {
      dd <- d
      for(s in unique(d$land_use)) {
        d <- dd[dd$land_use==s,]
        res[i, "start"] <- t
        res[i, "end"] <- t+(period-1)
        res[i,"land_use"] <- s
        if(!is.na(species_group)[1]) res[i,"spec_group"] <- "All"
        res[i, "N"] <- nrow(d)
        res[i, "MAR"] <- round(mean(d$SR), 2)
        try({
          res[i, "CIL"] <- round(t.test(d$SR)$conf.int[1], 2)
          res[i, "CIU"] <- round(t.test(d$SR)$conf.int[2], 2)
        })
        res[i, "REF"] <- round(mean(d$SR[d$SR >= quantile(d$SR, 0.75)]), 2)
        speciesID <- dat[dat$land_use==s & dat$pyear >= t & dat$pyear <= (t+(period-1)), "speciesID"]
        res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
        i <- i+1
      }			
    }  	
  }	
  
  ### Calculate MAR for the different species groups
  if(!is.na(species_group)[1]) {
    for(g in unique(species_group$group)) {
      
      for(t in start:(end-(period-1))) {
        d  <- surveys[surveys$pyear >= t & surveys$pyear <= (t+(period-1)), ]
        d$SR <- 0
        temp <- table(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g, "surveyID"])
        if(nrow(d) < length(temp)) message("The data.frame 'dat' contains observations from surveys that are not contained in 'surveys'.")
        d[names(temp),"SR"] <- temp
        d <- d[!is.na(d$coordID),]
        res[i, "start"] <- t
        res[i, "end"] <- t+(period-1)
        res[i,"land_use"] <- "All"
        res[i,"spec_group"] <- g
        res[i, "N"] <- nrow(d)
        res[i, "MAR"] <- round(mean(d$SR), 2)
        try({
          res[i, "CIL"] <- round(t.test(d$SR)$conf.int[1], 2)
          res[i, "CIU"] <- round(t.test(d$SR)$conf.int[2], 2)
        })
        res[i, "REF"] <- round(mean(d$SR[d$SR >= quantile(d$SR, 0.75)]), 2)
        speciesID <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g, "speciesID"]
        res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
        i <- i+1
        if(!is.na(land_use)[1]) {
          dd <- d
          for(s in unique(d$land_use)) {
            d <- dd[dd$land_use==s,]
            res[i, "start"] <- t
            res[i, "end"] <- t+(period-1)
            res[i,"land_use"] <- s
            res[i,"spec_group"] <- g
            res[i, "N"] <- nrow(d)
            res[i, "MAR"] <- round(mean(d$SR), 2)
            try({
              res[i, "CIL"] <- round(t.test(d$SR)$conf.int[1], 2)
              res[i, "CIU"] <- round(t.test(d$SR)$conf.int[2], 2)
            })
            res[i, "REF"] <- round(mean(d$SR[d$SR >= quantile(d$SR, 0.75)]), 2)
            speciesID <- dat[dat$land_use==s & dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g, "speciesID"]
            res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
            i <- i+1
          }			
        }  	
      }
      
    } # end of g-loop
  } # end of if()
  
  ### Order res and give it as output
  if(is.na(species_group)[1]) res <- res[order(res$land_use, res$start),]
  if(!is.na(species_group)[1]) res <- res[order(res$land_use, res$spec_group, res$start),]
  res
}
