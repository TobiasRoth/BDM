#--------------------------------------------------------------------------------------------------------------------
# Check whether one species is recorded several times during one survey
#--------------------------------------------------------------------------------------------------------------------
checkDouble <- 
function(dat) {
### Prepare and control data
	dat <- dat[,1:3]
	names(dat) <- c("coordID", "pyear", "speciesID")
	# Make unique ID for survey
	dat$surveyID <- paste(dat$coordID, "_", dat$pyear, sep="")
	# Remove not identified species
	dat <- dat[!is.na(dat$speciesID) & dat$speciesID != "", ]	
### Check whether one species is recorded several times during one survey
	temp <- apply(table(dat$surveyID, dat$speciesID)>1, 1, sum)
### Return results
	if(sum(temp)==0) res <- "No species is recorded two times during  a single survey."
	if(sum(temp)!=0) res <- names(temp[temp>0])
    res
}


#--------------------------------------------------------------------------------------------------------------------
# Check whether surveys are listed more than once in Kopfdaten
#--------------------------------------------------------------------------------------------------------------------
checkDoubleSurv <-
function(dat) {
  ### Drop all further columns except coordID and pyear
  dat<-dat[,1:2]
  ### Check for plots that were surveyed more than once within any five-year interval
  double.rows<-vector()
  for (i in 1:dim(dat)[1]) {
    same.coord<-which(dat[,1]==dat[i,1])
    if(length(same.coord)>1){
      differences<-abs(dat[same.coord[which(same.coord!=i)],2]-dat[i,2])
      if(min(differences)<5) c(double.rows,i)->double.rows
    }
  }
  ### Save rows containing plots that were surveyed more than once within five years
  temp <- dat[double.rows,]
  ### Check the number of surveys per year
  nobs <- data.frame(period=paste(min(dat[,2]):(max(dat[,2])-4), (min(dat[,2])+4):(max(dat[,2])), sep="-"))
  for (t in 1:nrow(nobs)) {
    ttt <- as.integer(substr(nobs[t,"period"], start=1, stop=4))
    nobs[t,"N"] <- sum(dat[,2] >= ttt & dat[,2] <= ttt+4)
  }
  ### Return results
  if(dim(temp)[1]==0) {
    cat("No double surveys in any five-year interval!", "\n\n", "Number of surveys per five-year intervall: \n")
    print(nobs)
  }
  if(dim(temp)[1]>0) {
    cat("Plots that were surveyed more than once within five years:","\n")
    print(temp[,1:2],row.names=F)
  }
}
