
#### Read Covid Data for Italy and Germany ####
#### Tobias Ruettenauer ####
#### 2020/ 03 / 16 ####

rm(list=ls())

### Load packages
library(rgdal)
library(spdep)
library(rgeos)
library(foreign)
library(GISTools)
library(cleangeo)
library(viridis)
library(magick)

library(extrafont)
loadfonts()


### Working Directory
setwd("C:/work/Forschung/Covid-19/02_Data")


### Italy data directory (https://github.com/pcm-dpc/COVID-19)
itad <- "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-province/dpc-covid19-ita-province.csv"


### Germany data directory (https://survstat.rki.de/) # Old data source
### Germany data directory (https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0) # New data source
ded <- "https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv"




###############################
### Prepare long data Italy ###
###############################

### Italian province data
italy.df <- read.table(itad,
                        header = TRUE, sep = ",", na.strings = "",
                        quote = "\"")
# Format date
italy.df$date <- as.Date(substr(italy.df$data, 1, 10), "%Y-%m-%d")

### Reorder
vars <- c("codice_provincia", "date", "denominazione_provincia", "sigla_provincia")
italy.df <- italy.df[c(vars, names(italy.df)[-which(names(italy.df) %in% vars)])]

italy.df <- italy.df[order(italy.df$codice_provincia, italy.df$date),]


### Total number of cases on 2020-20-14
sum(italy.df$totale_casi[italy.df$date == max(italy.df$date)])


### New cases per day per province
italy.df$new_cases <- ave(italy.df$totale_casi, 
                          by = italy.df$codice_provincia,
                          FUN = function(x) x - dplyr::lag(x))

sum(italy.df$new_cases[italy.df$date == max(italy.df$date)])

# Why negative "new cases"? Wrong test or moving between provinces?


### Rename ID
names(italy.df)[which(names(italy.df) == "codice_provincia")] <- "COD_PROV"
italy.df$COD_PROV <- as.character(italy.df$COD_PROV)



# #################################
# ### Prepare long data Germany ### ### Old data source
# #################################
# 
# germany_long.df <- reshape(germany.df[, -which(names(germany.df) %in% c("Anzahl", "Gesamt"))],
#                            varying = names(germany.df[, -c(1:3)]),
#                            direction = "long", idvar = "Kreis", sep = "_") 
# 
# names(germany_long.df)[which(names(germany_long.df) == "2020")] <- "daily_cases"
# rownames(germany_long.df) <- 1:nrow(germany_long.df)
# 
# 
# ### Add geo identifier
# # from https://gist.github.com/afternoon/1258046
# 
# id.df <- read.table("RKI_georef.txt",
#                     header = TRUE, sep = ",", na.strings = "",
#                     quote = "\'", skipNul = TRUE, strip.white = TRUE)
# names(id.df)[1] <- "id"
# 
# # Manually change
# germany_long.df$Kreis <- as.character(germany_long.df$Kreis)
# id.df$Kreis <- as.character(id.df$Kreis)
# 
# germany_long.df$Kreis <- gsub("\023 ", " ", germany_long.df$Kreis)
# 
# id.df$Kreis[id.df$Kreis == "LK Burgenland"] <- "LK Burgenlandkreis"
# id.df$Kreis[id.df$Kreis == "Heidekreis"] <- "LK Heidekreis"
# id.df$Kreis[id.df$Kreis == "LK Ludwigslust"] <- "LK Ludwigslust Parchim"
# id.df$Kreis[id.df$Kreis == "LK Mansfelder-Südharz"] <- "LK Mansfeld-Südharz"
# id.df$Kreis[id.df$Kreis == "LK Mecklenburg-Strelitz"] <- "LK Mecklenburgische Seenplatte"
# id.df$Kreis[id.df$Kreis == "LK Erftkreis"] <- "LK Rhein-Erft-Kreis"
# id.df$Kreis[id.df$Kreis == "LK Neuss"] <- "LK Rhein-Kreis Neuss"
# id.df$Kreis[id.df$Kreis == "LK Güstrow"] <- "LK Rostock"
# id.df$Kreis[id.df$Kreis == "LK Salzland"] <- "LK Salzlandkreis"
# id.df$Kreis[id.df$Kreis == "LK Ostvorpommern"] <- "LK Vorpommern Greifswald"
# id.df$Kreis[id.df$Kreis == "LK Nordvorpommern"] <- "LK Vorpommern Rügen"
# id.df$Kreis[id.df$Kreis == "SK Halle (Saale)"] <- "SK Halle"
# id.df$Kreis[id.df$Kreis == "SR Aachen"] <- "StädteRegion Aachen"
# 
# id.df$id[id.df$Kreis == "LK Ludwigslust Parchim"] <- 13076
# id.df$id[id.df$Kreis == "LK Mecklenburgische Seenplatte"] <- 13071
# id.df$id[id.df$Kreis == "LK Rostock"] <- 13072
# id.df$id[id.df$Kreis == "LK Vorpommern Greifswald"] <- 13075
# id.df$id[id.df$Kreis == "LK Vorpommern Rügen"] <- 13073
# 
# 
# # Merge
# 
# 
# germany_long.df <- merge(germany_long.df, id.df, by = "Kreis", all.x = TRUE)



#################################
### Prepare long data Germany ### 
#################################

### Germany daily RKI cases
germany_long.df <- read.table(ded, # enc: with encoding
                              header = TRUE, sep = ",", na.strings = "",
                              quote = "", skipNul = TRUE,
                              colClasses = c(IdLandkreis = "character"))
names(germany_long.df)[1] <- "id_bl"
names(germany_long.df)[which(names(germany_long.df) == "IdLandkreis")] <- "AGS"


### extract date
germany_long.df$date <- substr(germany_long.df$Meldedatum, 1, 10)
germany_long.df$date <- as.Date(germany_long.df$date, format = "%Y-%m-%d")

### extract age
germany_long.df$age_gr <- gsub("A", "", germany_long.df$Altersgruppe)
germany_long.df$age_gr <- gsub("-", "_", germany_long.df$age_gr)
germany_long.df$age_gr <- gsub("+", "", germany_long.df$age_gr, fixed = TRUE)



### set "unbekannt" NA
germany_long.df$age_gr[germany_long.df$age_gr == "unbekannt"] <- NA
germany_long.df$Geschlecht[germany_long.df$Geschlecht == "unbekannt"] <- NA


### total cases 
germany_long.df <- germany_long.df[order(germany_long.df$AGS, germany_long.df$date), ]

germany_long.df$daily_cases <- ave(germany_long.df$AnzahlFall,
                                   by = paste0(germany_long.df$AGS, germany_long.df$date),
                                   FUN = function(x) sum (x, na.rm = TRUE))
germany_long.df$daily_deaths <- ave(germany_long.df$AnzahlTodesfall,
                                   by = paste0(germany_long.df$AGS, germany_long.df$date),
                                   FUN = function(x) sum (x, na.rm = TRUE))


### reshape age group cases to wide
age.df <- germany_long.df[, c("AGS", "date", "age_gr", "AnzahlFall", "AnzahlTodesfall")]

# aggregate over gender
age.df$cases <- ave(germany_long.df$AnzahlFall,
                    by = paste0(germany_long.df$AGS, germany_long.df$date, germany_long.df$age_gr),
                    FUN = function(x) sum (x, na.rm = TRUE))
age.df$deaths <- ave(germany_long.df$AnzahlTodesfall,
                     by = paste0(germany_long.df$AGS, germany_long.df$date, germany_long.df$age_gr),
                     FUN = function(x) sum (x, na.rm = TRUE))
age.df[, c("AnzahlFall", "AnzahlTodesfall")] <- NULL
age.df <- unique(age.df)

age_w.df <- reshape(age.df, idvar = c("AGS", "date"), timevar = "age_gr",
                    direction = "wide", sep = "_age_")
age_w.df[is.na(age_w.df)] <- 0


### reshape gender cases to wide
sex.df <- germany_long.df[, c("AGS", "date", "Geschlecht", "AnzahlFall", "AnzahlTodesfall")]

# aggregate over age groups
sex.df$cases <- ave(germany_long.df$AnzahlFall,
                    by = paste0(germany_long.df$AGS, germany_long.df$date, germany_long.df$Geschlecht),
                    FUN = function(x) sum(x, na.rm = TRUE))
sex.df$deaths <- ave(germany_long.df$AnzahlTodesfall,
                     by = paste0(germany_long.df$AGS, germany_long.df$date, germany_long.df$Geschlecht),
                     FUN = function(x) sum(x, na.rm = TRUE))
sex.df[, c("AnzahlFall", "AnzahlTodesfall")] <- NULL
sex.df <- unique(sex.df)

sex_w.df <- reshape(sex.df, idvar = c("AGS", "date"), timevar = "Geschlecht",
                    direction = "wide", sep = "_sex_")
sex_w.df[is.na(sex_w.df)] <- 0



### Combine data
germany_long.df <- merge(germany_long.df, age_w.df, by = c("AGS", "date"),
                         all.x = TRUE, all.y = TRUE)
germany_long.df <- merge(germany_long.df, sex_w.df, by = c("AGS", "date"),
                         all.x = TRUE, all.y = TRUE)


### Drop report specific info and make unique
germany_long.df[, which(names(germany_long.df) %in% c("Geschlecht", "Altersgruppe", "age_gr", 
                                                      "ObjectId", "AnzahlFall", "AnzahlTodesfall",
                                                      "NeuerFall", "NeuerTodesfall"))] <- NULL

germany_long.df <- unique(germany_long.df)


### Cumulative cases by langkreis
germany_long.df <- germany_long.df[order(germany_long.df$AGS, germany_long.df$date), ]

germany_long.df$sum_cases <- ave(germany_long.df$daily_cases,
                                 by = germany_long.df$AGS,
                                 FUN = function(x) cumsum(x))
germany_long.df$sum_deaths <- ave(germany_long.df$daily_deaths,
                                 by = germany_long.df$AGS,
                                 FUN = function(x) cumsum(x))


### Cumulative cases by agegroup and langkreis
for(i in names(table(age.df$age_gr))){
  name1 <- paste0("cases_age_", i)
  name2 <- paste0("deaths_age_", i)
  newname1 <- paste0("sum_cases_", i)
  newname2 <- paste0("sum_deaths_", i)
  
  germany_long.df[, newname1] <- ave(germany_long.df[, name1],
                                   by = germany_long.df$AGS,
                                   FUN = function(x) cumsum(x))
  germany_long.df[, newname2] <- ave(germany_long.df[, name2],
                                    by = germany_long.df$AGS,
                                    FUN = function(x) cumsum(x))
}



############
### Save ###
############

### Save
save(italy.df, file = "Italy_covid.RData")
save(germany_long.df, file = "Germany_covid.RData")







