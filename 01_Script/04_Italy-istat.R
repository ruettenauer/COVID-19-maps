
#### Combine italy covid with istat ####
#### Tobias Ruettenauer ####
#### 2020/ 03 / 19 ####

rm(list=ls())

### Load packages
library(rgdal)
library(spatialreg)
library(spdep)
library(rgeos)
library(foreign)
library(readxl)
library(GISTools)
library(cleangeo)
library(viridis)
library(magick)
library(splm)
library(stringr)

library(extrafont)
loadfonts()


### Working Directory
setwd("C:/work/Forschung/Covid-19/02_Data")

### WHO Data directory (https://github.com/CSSEGISandData/COVID-19)
whod <- "C:/work/Forschung/Daten/COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

### Italy data directory (https://github.com/pcm-dpc/COVID-19)
itad <- "C:/work/Forschung/Daten/COVID-19-ita/dati-province/"


#################
### Load data ###
#################


load("Italy_covid.RData")


########################################
### Add italian poluation data  data ###
########################################
# Source: http://dati.istat.it/Index.aspx?QueryId=18460&lang=en


istat.df <- read.table("./Italian_population/DCIS_POPRES1_16032020110716829.csv", 
                       header = TRUE, sep = ",", na.strings = "",
                       quote = "\"")

names(istat.df)[1] <- "id"

### For now, keep only total cases (no cross-tab stats between indicators)
o1 <- as.integer(istat.df$Age == "total")
o2 <- as.integer(istat.df$Gender == "total")
o3 <- as.integer(istat.df$Marital.status == "total")

oo <- rowSums(cbind(o1, o2, o3))
oo <- which(oo >= 2)

istat.df <- istat.df[oo, ]


### Cut age into groups

istat.df$age_cat <- gsub(" .*$", "", istat.df$Age)
istat.df$age_cat[istat.df$age_cat == "total"] <- 999 
istat.df$age_cat <- as.numeric(istat.df$age_cat)

cuts <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 101, 999)

istat.df$age_cat <- cut(istat.df$age_cat, breaks = cuts, include.lowest = TRUE, right = FALSE)

### Aggregate by age
idvars <- c("id", "Territory")
vars <- c("SEXISTAT1", "Gender", "STATCIV2", "Marital.status", "age_cat")

istat.df <- istat.df[, c(idvars, vars, "Value")]

list <- as.list(istat.df[, c(idvars, vars)])

istat_agg.df <- aggregate(istat.df$Value, by = list, 
                          FUN = function(x) sum(x))


### Make wide

# Spread total
oo <- which(istat_agg.df$Gender == "total" & istat_agg.df$Marital.status == "total" & istat_agg.df$age_cat == "[101,999]")
istat_wide.df <- istat_agg.df[oo, c(idvars, "x")]
names(istat_wide.df)[which(names(istat_wide.df) == "x")] <- "total"

# Spread age
oo <- which(istat_agg.df$Gender == "total" & istat_agg.df$Marital.status == "total" & istat_agg.df$age_cat != "[101,999]")
tmp.df <- tidyr::spread(istat_agg.df[oo, c(idvars, "age_cat", "x")],
                        key = "age_cat", value = "x")
names(tmp.df)[-c(1:2)] <- paste0("age_", names(tmp.df)[-c(1:2)])  

istat_wide.df <- merge(istat_wide.df, tmp.df, by = idvars)

# Spread Marital status
oo <- which(istat_agg.df$Gender == "total" & istat_agg.df$Marital.status != "total" & istat_agg.df$age_cat == "[101,999]")
tmp.df <- tidyr::spread(istat_agg.df[oo, c(idvars, "Marital.status", "x")],
                        key = "Marital.status", value = "x")
names(tmp.df)[-c(1:2)] <- gsub(" ", "-", paste0("marstat_", names(tmp.df)[-c(1:2)]), fixed = TRUE)  
names(tmp.df)[-c(1:2)] <- gsub("/", "_", names(tmp.df)[-c(1:2)], fixed = TRUE)

istat_wide.df <- merge(istat_wide.df, tmp.df, by = idvars)

# Spread Gender
oo <- which(istat_agg.df$Gender != "total" & istat_agg.df$Marital.status == "total" & istat_agg.df$age_cat == "[101,999]")
tmp.df <- tidyr::spread(istat_agg.df[oo, c(idvars, "Gender", "x")],
                        key = "Gender", value = "x")
names(tmp.df)[-c(1:2)] <- gsub(" ", "-", paste0("sex_", names(tmp.df)[-c(1:2)]), fixed = TRUE)  

istat_wide.df <- merge(istat_wide.df, tmp.df, by = idvars)
  




##################
### Merge data ###
##################
# Use territory codes from http://wiki.scuola247.org/images/f/fe/Ripartizioni_regioni_province.txt

id.df <- read.table("Ripartizioni_regioni_province.txt", 
                    header = TRUE, sep = ";", na.strings = "-",
                    quote = "\"", nrow = 110, fill = TRUE, blank.lines.skip = TRUE)


### Reduce to relevant info

id.df <- id.df[, c("Codice.provincia", "Codice.NUTS3.2006", "Denominazione.provincia")]

names(id.df) <- c("COD_PROV", "id", "name")

id.df <- id.df[which(!is.na(id.df$COD_PROV)),]


### Manually correct 
id.df$id <- as.character(id.df$id)
id.df$name <- as.character(id.df$name)


id.df$id[which(id.df$name == "Monza e della Brianza")] <- "IT108"
id.df$id[which(id.df$name == "Fermo")] <- "IT109"
id.df$id[which(id.df$name == "Barletta-Andria-Trani")] <- "IT110"

# Add sud sardegna
id.df <- rbind(id.df, c(111, "IT111", "Sud Sardegna"))


### Merge IDs with istat data
istat_wide.df <- merge(istat_wide.df, id.df, by = "id", all.x = TRUE)

istat_wide.df <- istat_wide.df[which(!is.na(istat_wide.df$COD_PROV)), ]

### Merge with covid data

italy_pop.df <- merge(italy.df, istat_wide.df, by = "COD_PROV")


### Save
save(italy_pop.df, file = "Italy_covid_population.RData")



######################
### Load shapefile ###
######################

provinces.sp <- readOGR(dsn = "./Limiti01012019_g/ProvCM01012019_g",
                        layer = "ProvCM01012019_g_WGS84") 

provinces.sp <- spTransform(provinces.sp, CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))



italy.sp <- readOGR(dsn = "./Limiti01012019_g/RipGeo01012019_g",
                        layer = "RipGeo01012019_g_WGS84") 

italy.sp <- spTransform(italy.sp, CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


italy.sp <- gBuffer(italy.sp, byid = FALSE, width = 500)



#######################
### Prepare data ###
#######################


### Add are to data set

tmp <- data.frame(provinces.sp)[, c("COD_PROV", "Shape_Area")]
names(tmp) <- c("COD_PROV", "area")
tmp$area <- tmp$area / (1000 * 1000)

italy_pop.df <- merge(italy_pop.df, tmp, by = "COD_PROV")


### Share of age groups 

vars <- c("age_[0,10)", "age_[10,20)", "age_[20,30)", "age_[30,40)", "age_[40,50)", 
          "age_[50,60)", "age_[60,70)", "age_[70,80)", "age_[80,101)")

for(i in vars){
  j <- gsub("[", "per_", i, fixed = TRUE)
  j <- gsub(",", "_", j, fixed = TRUE)
  j <- gsub(")", "", j, fixed = TRUE)
  
  italy_pop.df[, j] <- italy_pop.df[, i] / italy_pop.df$total * 100
}


### Share of men

italy_pop.df$male_per <- italy_pop.df$sex_males / italy_pop.df$total * 100


# #######################
# ### listwise object ###
# #######################
# 
# it.nb <- poly2nb(provinces.sp, row.names = provinces.sp$COD_PROV)
# it.lw <- nb2listw(it.nb, style = "W")
# 
# 
# ### W multiplied over time
# W <- nb2mat(it.nb)
# WI <- W %x% diag(length(unique(italy_pop.df$date)))
# 
# WI2 <- W %x% diag((length(unique(italy_pop.df$date)) - 1))
# 
# 
# it_long.lw <- mat2listw(WI, style = "W")
# it_long2.lw <- mat2listw(WI2, style = "W")
# 
# 
# ############################################
# ### Spatial Autoregressive models ###
# ############################################
# 
# ### Drop all pending cases (not allocated to province)
# italy_red.df <- italy_pop.df[which(as.numeric(italy_pop.df$COD_PROV) < 900), ]
# 
# 
# ### Empty SAR model
# 
# sar0_pols.mod <- spml(totale_casi ~ 1, 
#                     data = italy_red.df, index = c("COD_PROV", "date"), listw = it.lw, 
#                     model = "pooling", effect = "individual",
#                     lag = TRUE, spatial.error = "none", tol.solve = 1e-10)
# summary(sar0_pols.mod)
# 
# 
# ### Pooled SAR
# sar1_pols.mod <- lagsarlm(totale_casi ~ total + area + male_per
#                           + age_per_0_10 + age_per_10_20 + age_per_30_40 + age_per_40_50 + 
#                             + age_per_50_60 + age_per_60_70 + age_per_70_80 + age_per_80_101, 
#                       data = italy_red.df, listw = it_long.lw,
#                       Durbin = FALSE, tol.solve = 1e-24 )
# summary(sar1_pols.mod)
# 
# impacts(sar1_pols.mod, listw = it_long.lw)
# 
# sar2_pols.mod <- lagsarlm(new_cases ~ total + area + male_per
#                           + age_per_0_10 + age_per_10_20 + age_per_30_40 + age_per_40_50 + 
#                             + age_per_50_60 + age_per_60_70 + age_per_70_80 + age_per_80_101,
#                           data = italy_red.df[italy_red.df$date != "2020-02-24", ], listw = it_long2.lw,
#                           Durbin = FALSE, tol.solve = 1e-24 )
# summary(sar2_pols.mod)
# 
# 
# ### Pooled SDM
# sdm1_pols.mod <- lagsarlm(totale_casi ~ total + area + male_per
#                           + age_per_0_10 + age_per_10_20 + age_per_30_40 + age_per_40_50 + 
#                             + age_per_50_60 + age_per_60_70 + age_per_70_80 + age_per_80_101, 
#                           data = italy_red.df, listw = it_long.lw,
#                           Durbin = TRUE, tol.solve = 1e-24 )
# summary(sdm1_pols.mod)





