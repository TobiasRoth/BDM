coordID2coord <- 
function(coordID) {
	data.frame(list(
		x = as.integer(substr(coordID, 1, 3))*1000,
		y = as.integer(substr(coordID, 4, 6))*1000
	))
}

#--------------------------------------------------------------------------------------------------------------------
# Get real species for unidentified species
#--------------------------------------------------------------------------------------------------------------------
simZA <- 
function(coordID, speciesID) {
  urne <- speciesID[!is.na(speciesID) & speciesID!=""]
  for(tt in unique(coordID)) {
    t.urne <- urne
		tspeciesID <- speciesID[!is.na(speciesID) & speciesID!="" & coordID==tt]
	  for(ttt in tspeciesID) t.urne <- t.urne[t.urne != ttt]
		n_ersetz <- sum((is.na(speciesID) | speciesID=="") & coordID==tt)
		arten_ersetz <- NULL
		if(n_ersetz>0){
		  for(jj in 1:n_ersetz) {
			  t.art <- sample(t.urne,1)
			  arten_ersetz <- c(arten_ersetz, t.art)
			  t.urne <- t.urne[t.urne!=t.art]
			} 
		  speciesID[(is.na(speciesID) | speciesID=="") & coordID==tt] <- arten_ersetz
		}
	}
  speciesID
}





