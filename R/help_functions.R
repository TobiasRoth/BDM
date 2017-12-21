#####################################################################################################################
# 'Summary' functions
#####################################################################################################################

#--------------------------------------------------------------------------------------------------------------------
# Function that gives some useful summary statistics of BDM data
#--------------------------------------------------------------------------------------------------------------------
summaryBDM <- 
function(dat) {
### Prepare and control data
	dat <- dat[,1:4]
	names(dat) <- c("coordID", "pyear", "speciesID", "speciesName")
	# Make unique ID for survey
	dat$surveyID <- paste(dat$coordID, "_", dat$pyear, sep="")
	
### How many surveys in total
	n_surveys <- length(unique(dat$surveyID))
		
### How many plots per year?
	temp <- table(dat$coordID, dat$pyear)
	n_plots <- apply(temp>0, 2, sum)
		
### Species richness per year
	SR_mean <- apply(temp, 2, function(x) mean(x[x>0]))
	SR_sd <- apply(temp, 2, function(x) sd(x[x>0]))
	
### How many species in total
	n_rec_species <- length(unique(dat$speciesID[!is.na(dat$speciesID) & !dat$speciesID==""]))
		
### How many unidentified species in total?
	temp <- (tapply(is.na(dat$speciesID) | dat$speciesID=="", dat$surveyID, sum) / table(dat$surveyID))*100
	unident_species  <- c(mean(temp), sd(temp))
	names(unident_species) <- c("mean", "sd")
	
### Frequency of species
	species <- data.frame(sort(table(dat$speciesName), decreasing=TRUE))
	names(species) <- c("n_records")

### Return results
	list(n_surveys=n_surveys, n_plots=n_plots, SR_mean =SR_mean, SR_sd =SR_sd, n_rec_species=n_rec_species, unident_species= unident_species, species=species)
}

# data(birdsLANAG)
# dat <- birdsLANAG[, c("coordID", "pyear", "speciesID", "speciesName")]
# MAR(dat)
#summaryBDM(dat)


#####################################################################################################################
# 'GIS' functions
#####################################################################################################################

#--------------------------------------------------------------------------------------------------------------------
# Convert UTM or Lon/Lat coordinates to Swiss grid
#--------------------------------------------------------------------------------------------------------------------
swissGrid <-
function(x, y = NA, UTM = FALSE, zone=NA, ...) {
	if(is.data.frame(x)) {
		y <- x[,2]
		x <- x[,1]
	}
	if(UTM) {
		tdat <- data.frame(list(X=x, Y=y))
		attr(tdat, "zone") <- zone
		attr(tdat, "projection") <- "UTM"
		temp <- PBSmapping::convUL(tdat, km=FALSE)
		x <- temp[,1]
		y <- temp[,2]
	}
	Hx <- (x*3600-26782.5)/10000
	Hy <- (y*3600-169028.66)/10000
	swissX <- 600072.37+211455.93*Hx-10938.51*Hx*Hy-0.36*Hx*Hy^2-44.54*Hx^3
	swissY <- 200147.07+308807.95*Hy+3745.25*Hx^2+76.63*Hy^2+119.79*Hy^3-194.56*Hx^2*Hy
	data.frame(list(x=swissX, y=swissY))
}

### UTM -> Swiss Grid
#swissGrid(393920, 5271690, projection="UTM", zone=32, km=FALSE)
### Lon/Lat -> Swiss Grid
#swissGrid(7.58901, 47.58989)


#####################################################################################################################
# 'Small' functions
#####################################################################################################################

#--------------------------------------------------------------------------------------------------------------------
# Probit an logit functions
#--------------------------------------------------------------------------------------------------------------------
probit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

#--------------------------------------------------------------------------------------------------------------------
# Transform plot Id to coordinates of Swiss grid
#--------------------------------------------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------------------------------------------
# Aggregates butterfly data from Swiss BDM
#--------------------------------------------------------------------------------------------------------------------
aggrButterflies <- function(dat) {
  print("Date of last modifications: 14.01.2015; based on 1180 BDM ArtenlisteMutterfileV9.xls")
  print("Important note: Wir haben eine Art (Erebia manto) neu aufgetrennt in die beiden Arten Erebia manto und Erebia bubastis! Das Auftrennen von Arten, von denen wir bereits Nachweise haben, sollte wenn möglich vermieden werden, denn sie sie löst allerlei Folgeprobleme und -arbeiten aus. Hier haben wir es trotzdem gemacht. Bis auf weiteres werden diese beiden Arten aber noch zusammengefasst.")

  dat[dat==29303] <- 10001	#Adscita statices-Komplex
  dat[dat==29304] <- 10001	#Adscita statices-Komplex
  dat[dat==29305] <- 10001	#Adscita statices-Komplex
  dat[dat==29307] <- 10001	#Adscita statices-Komplex
  dat[dat==29308] <- 10001	#Adscita statices-Komplex 
  
  dat[dat==31157] <- 31196 #Boloria napaea-Komplex
  dat[dat==31158] <- 31196 #Boloria napaea-Komplex
    
  dat[dat==29312] <- 10002 #Jordanita globularia-Komplex
  dat[dat==29313] <- 10002 #Jordanita globularia-Komplex
  dat[dat==29314] <- 10002 #Jordanita globularia-Komplex
  dat[dat==29315] <- 10002 #Jordanita globularia-Komplex
  
  dat[dat==31256] <- 31059 #Leptidea sinapis-Komplex
  
  dat[dat==31110] <- 31135 #Maculinea alcon-Komplex
  dat[dat==31113] <- 31135 #Maculinea alcon-Komplex
  
  dat[dat==31061] <- 31069 #Pieris napi-Komplex
  
  dat[dat==31062] <- 31070 #Pieris rapae-Komplex
	
  dat[dat==31010] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31011] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31012] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31013] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31014] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31015] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31017] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31020] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31021] <- 10004 #Pyrgus alveus-Komplex
  dat[dat==31022] <- 10004 #Pyrgus alveus-Komplex
	
	dat[dat==31018] <- 31027 #Pyrgus malvae-Komplex
	dat[dat==31019] <- 31027 #Pyrgus malvae-Komplex

  dat[dat==31263] <- 31223 #Erebia manto
  
	dat
	}




