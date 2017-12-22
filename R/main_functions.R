#--------------------------------------------------------------------------------------------------------------------
# MARdiff: Calculates differences in mean species richness between two periods
#--------------------------------------------------------------------------------------------------------------------
### Log of changes:
# - 03.11.11: Function implemented
MARdiff <-
function(dat, surveys, land_use=NA, species_group=NA, P1=c(2001,2005), P2=c(2006,2010) ) {
	
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
    land_use <- as.character(land_use)
    surveys$land_use <- land_use
	}

### Prepare 'species_group'
	if(!is.na(species_group)[1]) {
		if(!class(species_group)=="data.frame") stop("'species_group' should be a data.frame")
		if(ncol(species_group)!=2) stop("'species_group' should contain the two columns 'coordID' and 'group'")
		names(species_group) <- c("speciesID", "group")
		if(sum(table(species_group$speciesID)>1)>0) stop("Some species in 'species_group' are mentioned twice." )
    species_group$group <- as.character(species_group$group)
		row.names(species_group) <- as.character(species_group$speciesID)
		dat$group <- species_group[as.character(dat$speciesID),"group"]
	}
	
### Prepare 'P1' and 'P2'
	if(length(P1)!=2 | length(P2)!=2) stop("'P1' and 'P2' should give the starting and ending year of the two time periods to be analyced.")
	
### Prepare data.frame 'res' for the results and set counter 'i' to one
	res <- data.frame()
	i <- 1

### Calculate MAR for all records and every land use
	# Select the surveys of the first period
	surv  <- surveys[surveys$pyear >= P1[1] & surveys$pyear <= P1[2], ]
	if(sum(table(surv$coordID)>1)>0) stop("During the first period a plot was surveyd two times.")
	row.names(surv) <- surv$coordID
	temp <- surveys[surveys$pyear >= P2[1] & surveys$pyear <= P2[2], ]
	if(sum(table(temp$coordID)>1)>0) stop("During the second period a plot was surveyd two times.")
	surv[as.character(temp$coordID), "pyear2"] <- temp$pyear 
	if(!is.na(land_use)[1]) surv[as.character(temp$coordID), "land_use2"] <- temp$land_use
	temp <- sum(is.na(surv$pyear) | is.na(surv$pyear2))
	if(temp>0) message(paste("For", temp, "plot(s) only one survey was available and therefore dumped."))
	surv <- surv[!is.na(surv$pyear) & !is.na(surv$pyear2),]
	surv$surveyID2 <- paste(surv$coordID, "_", surv$pyear2, sep="")
	# Get changes in land_use
	if(!is.na(land_use)[1]) land_use_change <- table(surv$land_use, surv$land_use2) 
	if(is.na(land_use)[1]) land_use_change <- "land use not provided"
	# Get species richness from first Period
	temp <- table(dat$surveyID)
	surv$SR1 <- temp[surv$surveyID]
	surv$SR2 <- temp[surv$surveyID2]
	surv$SR1[is.na(surv$SR1)] <- 0
	surv$SR2[is.na(surv$SR2)] <- 0
	res[i,"land_use"] <- "all"
	if(!is.na(species_group)[1]) res[i,"spec_group"] <- "all"
	res[i, "N"] <- nrow(surv)
	res[i, "SR1"] <- round(mean(surv$SR1),2)
	res[i, "SR2"] <- round(mean(surv$SR2),2)
	res[i, "SR_Dif"] <- round(mean(surv$SR2-surv$SR1),2)
	try({
		res[i, "CIL"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[1],2)
		res[i, "CIU"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[2],2)
		temp <- t.test(surv$SR2-surv$SR1)$p.value
		res[i, "P.value"] <- ifelse(temp<0.001, "<0.001", round(temp ,3))
	})
	i <- i+1
	if(!is.na(land_use)[1]) {
	tsurv <- surv
	for(s in unique(land_use)) {
		surv <- tsurv[tsurv$land_use==s & tsurv$land_use2==s, ]
		res[i,"land_use"] <- s
		if(!is.na(species_group)[1]) res[i,"spec_group"] <- "all"
		res[i, "N"] <- nrow(surv)
		res[i, "SR1"] <- round(mean(surv$SR1),2)
		res[i, "SR2"] <- round(mean(surv$SR2),2)
		res[i, "SR_Dif"] <- round(mean(surv$SR2-surv$SR1),2)
		try({
			res[i, "CIL"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[1],2)
			res[i, "CIU"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[2],2)
			temp <- t.test(surv$SR2-surv$SR1)$p.value
			res[i, "P.value"] <- ifelse(temp<0.001, "<0.001", round(temp ,3))
		})
		i <- i+1
		}			
	}  	
	
### Calculate MARdiff for the different species groups
	if(!is.na(species_group)[1]) {
	for(g in unique(species_group$group)) {
		# Select the surveys of the first period
		surv  <- tsurv
		surv$SR1 <- 0
		surv$SR2 <- 0
		temp <- table(dat[species_group[as.character(dat$speciesID), "group"]==g & !is.na(dat$speciesID) & dat$speciesID!="", "surveyID"])
		surv$SR1 <- temp[surv$surveyID]
		surv$SR2 <- temp[surv$surveyID2]
		surv$SR1[is.na(surv$SR1)] <- 0
		surv$SR2[is.na(surv$SR2)] <- 0
		res[i,"land_use"] <- "all"
		res[i,"spec_group"] <- g
		res[i, "N"] <- nrow(surv)
		res[i, "SR1"] <- round(mean(surv$SR1),2)
		res[i, "SR2"] <- round(mean(surv$SR2),2)
		res[i, "SR_Dif"] <- round(mean(surv$SR2-surv$SR1),2)
		try({
			res[i, "CIL"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[1],2)
			res[i, "CIU"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[2],2)
			temp <- t.test(surv$SR2-surv$SR1)$p.value
			res[i, "P.value"] <- ifelse(temp<0.001, "<0.001", round(temp ,3))
		})
		i <- i+1
		if(!is.na(land_use)[1]) {
		tsurv <- surv
		for(s in unique(land_use)) {
			surv <- tsurv[tsurv$land_use==s & tsurv$land_use2==s, ]
			res[i,"land_use"] <- s
			res[i,"spec_group"] <- g
			res[i, "N"] <- nrow(surv)
			res[i, "SR1"] <- round(mean(surv$SR1),2)
			res[i, "SR2"] <- round(mean(surv$SR2),2)
			res[i, "SR_Dif"] <- round(mean(surv$SR2-surv$SR1),2)
			try({
				res[i, "CIL"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[1],2)
				res[i, "CIU"] <- round(t.test(surv$SR2-surv$SR1)$conf.int[2],2)
				temp <- t.test(surv$SR2-surv$SR1)$p.value
				res[i, "P.value"] <- ifelse(temp<0.001, "<0.001", round(temp ,3))
			})
			i <- i+1
		}			
		}  		
	} # end of g-loop
	} # end of if()

### Order res and give it as output
	if(is.na(species_group)[1]) {res <- res[order(res$land_use),]}
	if(!is.na(species_group)[1]) res <- res[order(res$land_use, res$spec_group),]
	list(res=res, land_use_change=land_use_change)
	
}
# data(birdsLANAG)
# dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
# surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
# surveys$pyear <- substr(surveys$coordID, 8,12)
# surveys$coordID <- substr(surveys$coordID, 1,6)
# MARdiff(dat, surveys)

#--------------------------------------------------------------------------------------------------------------------
# Z12: Calculates species disimilarity between plots from a given stratum.
#--------------------------------------------------------------------------------------------------------------------
### Log of changes:
# - 03.11.2011: Function implemented
# - 16.11.2011: Implemented simulation to account for not identified species
# - 17.11.2011: Possible indicator values are between 0 and 100 and not 0 and 1.
Z12 <-
function(dat, surveys, land_use=NA, species_group=NA, start=2001, end=2010, N_sim=0, signif=2, period=5, method="simpson") {
	
### Select relevant data from 'dat' and make unique ID for survey	
	if(class(dat) != "data.frame") stop("'dat' should be a data.frame")
	if(ncol(dat)<3) stop("'dat' should contain three columns")
	dat <- dat[,1:3]
	names(dat) <- c("coordID", "pyear", "speciesID")
	dat$surveyID <- paste(dat$coordID, "_", dat$pyear, sep="")
	temp <- sum(is.na(dat$speciesID) | dat$speciesID == "")
	if(temp > 0 & N_sim==0) message(paste("The", temp, "records of unidentified species are droped from 'dat'." ))
	if(N_sim == 0) dat <- dat[!(is.na(dat$speciesID) | dat$speciesID == ""),]
	dat$oc <- 1
	
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

### Calculate disimilarity index for all records and every land use
	for(t in start:(end-(period-1))) {
		d  <- surveys[surveys$pyear >= t & surveys$pyear <= (t+(period-1)), ]
		d$SR <- 0
		temp <- table(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)), "surveyID"])
		if(nrow(d) < length(temp)) stop("The data.frame 'dat' contains observations from surveys that are not contained in 'surveys'.")
		d[names(temp),"SR"] <- temp
		d <- d[!is.na(d$coordID),]
		res[i, "start"] <- t
		res[i, "end"] <- t+(period-1)
		res[i,"land_use"] <- "All"
		if(!is.na(species_group)[1]) res[i,"spec_group"] <- "All"
		res[i, "N"] <- nrow(d)
		res[i, "N0"] <- sum(d$SR==0)
		if(N_sim == 0) {
      try({
      ttt <- simba::sim(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)), c("coordID", "speciesID","oc")], method=method, listin=T, listout=T)
		  res[i, "Z12"] <- round(mean(ttt[,method])*100, signif)
			res[i, "CIL"] <- round(t.test(ttt[,method])$conf.int[1]*100, signif)
			res[i, "CIU"] <- round(t.test(ttt[,method])$conf.int[2]*100, signif) })
		}
  	if(N_sim > 0) {
      tres <- numeric(N_sim)
      try({
        for(l in 1:N_sim) {
         simdat <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)), c("coordID", "speciesID","oc")]
          simdat$speciesID <- simZA(simdat$coordID, simdat$speciesID)
          ttt <- simba::sim(simdat, method=method, listin=T, listout=T)
          tres[l] <- mean(ttt[,method])
        }
	  	  res[i, "Z12"] <- round(mean(tres)*100, signif)
  		  res[i, "CIL"] <- round(min(tres)*100, signif)
			  res[i, "CIU"] <- round(max(tres)*100, signif)
      })
  	}
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
				res[i, "N0"] <- sum(d$SR==0)
        if(N_sim == 0) {
          try({
				    ttt <- simba::sim(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & surveys[as.character(dat$surveyID), "land_use"]==s, c("coordID", "speciesID","oc")], method=method, listin=T, listout=T)
  	        res[i, "Z12"] <- round(mean(ttt[,method])*100, signif)
  		      res[i, "CIL"] <- round(t.test(ttt[,method])$conf.int[1]*100, signif)
			      res[i, "CIU"] <- round(t.test(ttt[,method])$conf.int[2]*100, signif) })
				}
        if(N_sim > 0) {
          tres <- numeric(N_sim)
          try({
            for(l in 1:N_sim) {
              simdat <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & surveys[as.character(dat$surveyID), "land_use"]==s, c("coordID", "speciesID","oc")]
              simdat$speciesID <- simZA(simdat$coordID, simdat$speciesID)
              ttt <- simba::sim(simdat, method=method, listin=T, listout=T)
              tres[l] <- mean(ttt[,method])
            }
    	      res[i, "Z12"] <- round(mean(tres)*100, signif)
  		      res[i, "CIL"] <- round(min(tres)*100, signif)
			      res[i, "CIU"] <- round(max(tres)*100, signif)
          })
		    }
        speciesID <- dat[dat$land_use==s & dat$pyear >= t & dat$pyear <= (t+(period-1)), "speciesID"]
				res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
        i <- i+1
			} # end land_use-loop			
		} # end if for land_use calculation
	} # end period-loop 

### Calculate Z12 for the different species groups
	if(!is.na(species_group)[1]) {
	for(g in unique(species_group$group)) {
			
	for(t in start:(end-(period-1))) {
		d  <- surveys[surveys$pyear >= t & surveys$pyear <= (t+(period-1)), ]
		d$SR <- 0
		temp <- table(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g, "surveyID"])		
		if(nrow(d) < length(temp)) stop("The data.frame 'dat' contains observations from surveys that are not contained in 'surveys'.")
		d[names(temp),"SR"] <- temp
		d <- d[!is.na(d$coordID),]
		res[i, "start"] <- t
		res[i, "end"] <- t+(period-1)
		res[i,"land_use"] <- "All"
		res[i,"spec_group"] <- g
		res[i, "N"] <- nrow(d)
		res[i, "N0"] <- sum(d$SR==0)
    if(N_sim == 0) {
      try({
		    ttt <- simba::sim(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g & !is.na(dat$speciesID) & dat$speciesID !="", c("coordID", "speciesID","oc")], method=method, listin=T, listout=T)
  	    res[i, "Z12"] <- round(mean(ttt[,method])*100, signif)
    		res[i, "CIL"] <- round(t.test(ttt[,method])$conf.int[1]*100, signif)
	  		res[i, "CIU"] <- round(t.test(ttt[,method])$conf.int[2]*100, signif) })
		}
    if(N_sim > 0) {
      tres <- numeric(N_sim)
      try({
        for(l in 1:N_sim) {
          simdat <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g & !is.na(dat$speciesID) & dat$speciesID !="", c("coordID", "speciesID","oc")]
          simdat$speciesID <- simZA(simdat$coordID, simdat$speciesID)
          ttt <- simba::sim(simdat, method=method, listin=T, listout=T)
  	      tres[l] <- mean(ttt[,method])
        }
    	  res[i, "Z12"] <- round(mean(tres)*100, signif)
  		  res[i, "CIL"] <- round(min(tres)*100, signif)
			  res[i, "CIU"] <- round(max(tres)*100, signif)
      }, silent=TRUE)
		}
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
				res[i,"spec_group"] <- g
				res[i, "N"] <- nrow(d)
				res[i, "N0"] <- sum(d$SR==0)
        if(N_sim == 0) {
          try({
            ttt <- simba::sim(dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g & !is.na(dat$speciesID) & dat$speciesID !="", c("coordID", "speciesID","oc")], method=method, listin=T, listout=T)
    	      res[i, "Z12"] <- round(mean(ttt[,method])*100, signif)
  		      res[i, "CIL"] <- round(t.test(ttt[,method])$conf.int[1]*100, signif)
			      res[i, "CIU"] <- round(t.test(ttt[,method])$conf.int[2]*100, signif) })
				}
        if(N_sim > 0) {
          tres <- numeric(N_sim)
          try({
            for(l in 1:N_sim) {
              simdat <- dat[dat$pyear >= t & dat$pyear <= (t+(period-1)) & dat$group==g & !is.na(dat$speciesID) & dat$speciesID !="", c("coordID", "speciesID","oc")]
              simdat$speciesID <- simZA(simdat$coordID, simdat$speciesID)
              ttt <- simba::sim(simdat, method=method, listin=T, listout=T)
              tres[l] <- mean(ttt[,method])
            }
    	      res[i, "Z12"] <- round(mean(tres)*100, signif)
  		      res[i, "CIL"] <- round(min(tres)*100, signif)
			      res[i, "CIU"] <- round(max(tres)*100, signif)
          })
  		  }			
  			speciesID <- dat[dat$land_use==s & dat$pyear >= t & dat$pyear <= (t+(period-1)), "speciesID"]
				res[i,"N_species"] <- length(unique(speciesID[!is.na(speciesID)]))
        i <- i+1
		  } # end species-group-loop  	
	  }
  
	} # end of t-loop
	} # end of g-loop
	} # end of if()

### Order res and give it as output
	if(is.na(species_group)[1]) res <- res[order(res$land_use, res$start),]
	if(!is.na(species_group)[1]) res <- res[order(res$land_use, res$spec_group, res$start),]
	res
	
}

# data(birdsLANAG)
# dat <- birdsLANAG[, c("coordID", "pyear", "speciesID")]
# surveys <- data.frame(list(coordID=unique(paste(dat$coordID, "_", dat$pyear, sep=""))))
# surveys$pyear <- substr(surveys$coordID, 8,12)
# surveys$coordID <- substr(surveys$coordID, 1,6)
# Z12(dat, surveys, start=2001, end=2010)

#--------------------------------------------------------------------------------------------------------------------
# KI: Calculates trend index for several species groups
#--------------------------------------------------------------------------------------------------------------------
### Log of changes:
# - 11.11.11: Function implemented
KI <-
function(dat,  base, signif=0, Z12=TRUE) {
### Test 'dat' and ,'base'
  if(length(dat) != length(base)) stop("'dat' and 'base' should be of the same length.")
  if(sum(names(dat) != names(base))>0) stop("The names of 'dat' and 'base' should correspond.")  
  res <- dat[[1]][,1:(which(names(dat[[1]])=="N")-1)]
  name_sg <- names(dat)
  number_sg <- length(name_sg)  
  for(i in 1:number_sg) {
    if(!Z12) res[,paste("MAR", i, sep="")] <- dat[[i]]$MAR/base[i]
    if(Z12) res[,paste("Z12", i, sep="")] <- dat[[i]]$Z12/base[i]
  }
  if(!Z12) { res$KI <- round(
  apply(res[,paste("MAR", 1:number_sg, sep="")], 1, 
        function(x) {
          prodx <- 1
          for(t in 1:length(x)) prodx <-  prodx * x[t]
          (prodx^(1/length(x))) * 100
        }
  ), signif) }
  if(Z12) { res$KI <- round(
  apply(res[,paste("Z12", 1:number_sg, sep="")], 1, 
        function(x) {
          prodx <- 1
          for(t in 1:length(x)) prodx <-  prodx * x[t]
          (prodx^(1/length(x))) * 100
        }
  ), signif) }
  res  
}
# data(LANAG)
# sur <- LANAG$sur
# surBU <- LANAG$surBu
# bi <- LANAG$bi
# pl <- LANAG$pl
# mo <- LANAG$mo
# bu <- LANAG$bu
# resmar <- list(
#   bi=MAR(bi, sur[!is.na(sur$yearBi),1:2], sur[!is.na(sur$yearBi),"land_use_bi"], start=1996, end=2011),
#   pl=MAR(pl, sur[!is.na(sur$yearPf),1:2], sur[!is.na(sur$yearPf),"land_use"], start=1996, end=2011),
#   mo=MAR(mo, sur[!is.na(sur$yearMo),1:2], sur[!is.na(sur$yearMo),"land_use"], start=1996, end=2011),
#   bu=MAR(bu, surBU[!is.na(surBU$yearBu),1:2], surBU[!is.na(surBU$yearBu),"land_use"], start=1996, end=2011)
# )
# resmar$bi <- resmar$bi[resmar$bi$land_use!="keine_Hauptnutzung",]
# resmar$bu <- resmar$bu[resmar$bu$land_use!="keine_Hauptnutzung",]
# base <- c(11.7, 14.0, 6.6, 6.1)
# names(base) <- c("bi", "pl", "mo", "bu")
# KI(resmar, base)

#--------------------------------------------------------------------------------------------------------------------
# Z8: Calculates population trends of single species
#--------------------------------------------------------------------------------------------------------------------
### Log of changes:
# - 11.11.2011: Function implemented
# - 12.11.2011: Poisson GLMM for trend in abundance of occupied sites
# - 18.11.2011: McNemar test added
# - 19.11.2011: McNemar test for 5 year period and 10 year period
# - 04.03.2015: Intercept for GLM fixed to year 2010
# - 04.03.2015: Signifikantstests wird nun auch für sehr seltene Arten durchgeführt. Fehlermeldungen werden aber abgefangen.
# - 04.03.2015: Trendberechnung ohne Quadrat in GLM
Z8 <-
function(dat, surveys, occ=TRUE, abundance=FALSE, McNemar5=FALSE, McNemar10=FALSE, delNull = FALSE) {

### Prepare 'dat' 
  if(class(dat) != "data.frame") stop("'dat' should be a data.frame")
  if((occ | McNemar5 | McNemar10) & ncol(dat)<4) stop("'dat' should contain at least four columns")
  if(abundance & ncol(dat)<5) stop("If trend from abundance shoud be calculated then 'dat' should contain at least five columns")
	if(ncol(dat)==4) {
	  names(dat) <- c("coordID", "pyear", "speciesID", "speciesName")
	}
  if(ncol(dat)>4) {
	  names(dat) <- c("coordID", "pyear", "speciesID", "speciesName", "Ind")
	}
	dat$surveyID <- paste(dat$coordID, "_", dat$pyear, sep="")
  temp <- sum(is.na(dat$speciesID) | dat$speciesID == "")
	if(temp > 0) message(paste("The", temp, "records of unidentified species are droped from 'dat'." ))
	dat <- dat[!(is.na(dat$speciesID) | dat$speciesID == ""),]

### Prepare 'surveys'
	if(class(surveys) != "data.frame") stop("'surveys' should be a data.frame")
	if(ncol(surveys)!=3) stop("'surveys' should contain the tree columns 'coordID',  'pyear' and 'year'")
	names(surveys) <- c("coordID", "pyear", "year")
	surveys$surveyID <- paste(surveys$coordID, "_", surveys$pyear, sep="")
	if(sum(table(surveys$surveyID)>1)>0) stop("There are more than one entry for a survey. Check in 'survey'")
	row.names(surveys) <- surveys$surveyID
  surveys$y <- (surveys$year-2010)/10
  tt <- sum(is.na(match(dat$surveyID, surveys$surveyID)))
  if(tt>0) message(paste("There are observations in dat that are not from the provided surveys. These", tt, "records are removed."))
  dat <- dat[!is.na(match(dat$surveyID, surveys$surveyID)),]

### Prepare data.frame 'res' for the results
	res <- data.frame(list(speciesID=sort(unique(dat$speciesID))))
  row.names(res) <- res$speciesID
  temp <- tapply(as.character(dat$speciesName), dat$speciesID, function(x) x[1])
	res[names(temp),"speciesName"] <- temp

### Prepare data.frame with sites
  sites <- data.frame(list(coordID=unique(surveys$coordID)))
  row.names(sites) <- sites$coordID
  endyear <- max(surveys$pyear)
  sites$OK1 <- FALSE
  sites[as.character(surveys[surveys$pyear>=endyear-14 & surveys$pyear<=endyear-10, "coordID"]),"OK1"] <- TRUE
  sites$OK2 <- FALSE
  sites[as.character(surveys[surveys$pyear>=endyear-9 & surveys$pyear<=endyear-5, "coordID"]),"OK2"] <- TRUE
  sites$OK3 <- FALSE
  sites[as.character(surveys[surveys$pyear>=endyear-4 & surveys$pyear<=endyear, "coordID"]),"OK3"] <- TRUE
  N <- sum(sites$OK2 & sites$OK3)

###  Start loop to calculate trend for each species
  for(i in 1:nrow(res)) {
    sites$Occ1 <- 0
    sites[as.character(dat[dat$pyear>=endyear-9 & dat$pyear<=endyear-5 & dat$speciesID==res$speciesID[i], "coordID"]),"Occ1"] <- 1        
    sites$Occ2 <- 0
    sites[as.character(dat[dat$pyear>=endyear-4 & dat$pyear<=endyear & dat$speciesID==res$speciesID[i], "coordID"]),"Occ2"] <- 1            
    res[i, "Occ_1"] <- round(100*sum(sites$OK2 & sites$OK3 & sites$Occ1 == 1)/N, 1)
    res[i, "Occ_2"] <- round(100*sum(sites$OK2 & sites$OK3 & sites$Occ2 == 1)/N, 1)
    res[i, "N_const"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1 == 1 & sites$Occ2 == 1)
    res[i, "N_new"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1 == 0 & sites$Occ2 == 1)
    res[i, "N_lost"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1 == 1 & sites$Occ2 == 0)
    surveys$Occ <- 0
    surveys[dat[dat$speciesID==res$speciesID[i], "surveyID"], "Occ"] <- 1 
    if(sum(names(dat) == "Ind", na.rm = TRUE) == 1) {
      surveys$Ind <- 0
      surveys[dat[dat$speciesID==res$speciesID[i], "surveyID"], "Ind"] <- dat[dat$speciesID==res$speciesID[i], "Ind"]
      sites$Ind1 <- 0
      sites[as.character(dat[dat$pyear>=endyear-9 & dat$pyear<=endyear-5 & dat$speciesID==res$speciesID[i], "coordID"]),"Ind1"] <- dat[dat$pyear>=endyear-9 & dat$pyear<=endyear-5 & dat$speciesID==res$speciesID[i], "Ind"]        
      sites$Ind2 <- 0
      sites[as.character(dat[dat$pyear>=endyear-4 & dat$pyear<=endyear & dat$speciesID==res$speciesID[i], "coordID"]),"Ind2"] <- dat[dat$pyear>=endyear-4 & dat$pyear<=endyear & dat$speciesID==res$speciesID[i], "Ind"]           
      res[i, "Ind_1"] <- round(mean(sites[sites$OK2 & sites$OK3 & (sites$Ind1 > 0 | sites$Ind2 > 0), "Ind1"]),2)
      res[i, "Ind_2"] <- round(mean(sites[sites$OK2 & sites$OK3 & (sites$Ind1 > 0 | sites$Ind2 > 0), "Ind2"]),2)      
    }
    ### Calculate trend from presence absence data    
    if(occ) {try({
      mod <- glmer(Occ ~ y + (1|coordID), data=surveys, family="binomial")
      res[i, "B_GLMM_L"] <- round(summary(mod)$coefficients[2,1], 3)
      res[i, "B_GLMM_L_p"] <- round(summary(mod)$coefficients[2,4], 3) #ifelse(summary(mod)$coefficients[2,4]<0.001, "<0.001", round(summary(mod)$coefficients[2,4], 3))
    })}
    
    ### Calcualte trend from abundance data      
    if(abundance) {try({
      surveys$Ind <- 0
      surveys[dat[dat$speciesID==res$speciesID[i], "surveyID"], "Ind"] <- dat[dat$speciesID==res$speciesID[i], "Ind"]
      if(delNull) {  
        sel <- tapply(surveys$Ind, surveys$coordID, sum)
        sel <- !is.na(match(surveys$coordID, names(sel)))
        mod <- glmer(Ind ~ y + (1|y) + (1|coordID), data=surveys[sel,], family="poisson")
      }
      if(!delNull) {
        mod <- glmer(Ind ~ y + (1|y) + (1|coordID), data=surveys, family="poisson")
      }
      res[i, "P_GLMM_L"] <- round(summary(mod)$coefficients[2,1], 3)
      res[i, "P_GLMM_L_p"] <- round(summary(mod)$coefficients[2,4], 3) #ifelse(summary(mod)$coefficients[2,4]<0.001, "<0.001", round(summary(mod)$coefficients[2,4], 3))
    })}
    
    ### Calculate trend using McNemar-test      
    if(McNemar5) {try({
      sites$Occ1 <- 0
      sites[as.character(dat[dat$pyear>=endyear-9 & dat$pyear<=endyear-5 & dat$speciesID==res$speciesID[i], "coordID"]),"Occ1"] <- 1        
      sites$Occ2 <- 0
      sites[as.character(dat[dat$pyear>=endyear-4 & dat$pyear<=endyear & dat$speciesID==res$speciesID[i], "coordID"]),"Occ2"] <- 1        
      res[i, "N5"] <- sum(sites$OK2 & sites$OK3)
      res[i, "Occ5"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1==1)
      res[i, "Increase5"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1==0 & sites$Occ2==1)
      res[i, "Decrease5"] <- sum(sites$OK2 & sites$OK3 & sites$Occ1==1 & sites$Occ2==0)
      res[i, "McNemar_p5"] <- p_wert <- ((abs(res[i, "Decrease5"]-res[i, "Increase5"])-1)^2)/(res[i, "Increase5"]+res[i, "Decrease5"])
      if(p_wert!=Inf) res[i, "McNemar_p5"] <- 1- pchisq(p_wert, 1) #ifelse(1-pchisq(p_wert, 1)<0.001, "<0.001", round(1-pchisq(p_wert, 1),3))    
    })}
    if(McNemar10) {try({
      sites$Occ1 <- 0
      sites[as.character(dat[dat$pyear>=endyear-14 & dat$pyear<=endyear-10 & dat$speciesID==res$speciesID[i], "coordID"]),"Occ1"] <- 1        
      sites$Occ2 <- 0
      sites[as.character(dat[dat$pyear>=endyear-4 & dat$pyear<=endyear & dat$speciesID==res$speciesID[i], "coordID"]),"Occ2"] <- 1        
      res[i, "N10"] <- sum(sites$OK1 & sites$OK3)
      res[i, "Occ10"] <- sum(sites$OK1 & sites$OK3 & sites$Occ1==1)
      res[i, "Increase10"] <- sum(sites$OK1 & sites$OK3 & sites$Occ1==0 & sites$Occ2==1)
      res[i, "Decrease10"] <- sum(sites$OK1 & sites$OK3 & sites$Occ1==1 & sites$Occ2==0)
      p_wert <- ((abs(res[i, "Decrease10"]-res[i, "Increase10"])-1)^2)/(res[i, "Increase10"]+res[i, "Decrease10"])
      if(p_wert!=Inf) res[i, "McNemar_p10"] <- 1- pchisq(p_wert, 1) #ifelse(1-pchisq(p_wert, 1)<0.001, "<0.001", round(1-pchisq(p_wert, 1),3))    
    })}      
  }

### Order res and give it as output
	list(res=res, N=N)	
}


