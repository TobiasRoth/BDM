#------------------------------------------------------------------------------------------------------
# Example data for BDM package
#------------------------------------------------------------------------------------------------------
rm(list = ls(all = TRUE))

# Libraries
library(tidyverse)

# Connection to data base
db <- src_sqlite(path = "~/Documents/Dropbox/DB-HW.db", create = FALSE)

#------------------------------------------------------------------------------------------------------
# Butterfly surveys from one site
#------------------------------------------------------------------------------------------------------
BDMZ7TF <- inner_join(tbl(db, "KD_Z7") %>%  filter(aID_STAO == 599222) %>% select(aID_KD, aID_STAO, yearBu),
                      tbl(db, "TF") %>% filter(!is.na(aID_SP)) %>% select(aID_KD, aID_SP, Ind)) %>% 
  left_join(tbl(db, "arten") %>% filter(TF) %>% transmute(aID_SP = aID_SP, genus = Gattung, species = Art)) %>% 
  transmute(siteID = aID_STAO + 102,
            year = yearBu,
            speciesID = aID_SP,
            genus = genus,
            species = species,
            number = Ind) %>% 
  data.frame

devtools::use_data(BDMZ7TF, overwrite = TRUE)
