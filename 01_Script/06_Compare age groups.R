
#### Combine RKI covid with INKAR ####
#### Tobias Ruettenauer ####
#### 2020/ 03 / 22 ####

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




#################
### Load data ###
#################


load("Germany_covid_full.RData")


### Combine Berlin to one County
berid <- unique(germany_covid.df$AGS[which(germany_covid.df$BL == "Berlin")])

ber.df <- germany_covid.df[germany_covid.df$AGS %in% berid, ]

drop <- c("AGS", "date", "SDV_RS", "GEN", "NUTS", "BL", "BL_ID", "Landkreis", "Meldedatum", "Datenstand")
ber.df$AGS <- "11000"

ber_agr.df <- aggregate(ber.df[, -which(names(ber.df) %in% drop)], 
                        by = list(AGS = ber.df$AGS, date = ber.df$date),
                        FUN = function(x) sum(x))

# Re-merge
germany_covid_agr.df <- dplyr::bind_rows(germany_covid.df[-which(germany_covid.df$AGS %in% berid), ], ber_agr.df)





######################
### Load shapefile ###
######################
# using RKI shapefile (https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/917fc37a709542548cc3be077a786c17_0)

germany.sp <- readOGR(dsn = "./RKI_Corona_Landkreise",
                      layer = "RKI_Corona_Landkreise") 
proj4string(germany.sp)
germany.sp <- spTransform(germany.sp, CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

germany.sp$AGS <- as.character(germany.sp$AGS)
germany.sp$RS <- as.character(germany.sp$RS)
germany.sp$AGS[which(is.na(germany.sp$AGS))] <- germany.sp$RS[which(is.na(germany.sp$AGS))]

### Aggregate Berlin

germany_agr.sp <- germany.sp
germany_agr.sp$AGS <- as.character(germany_agr.sp$AGS)
germany_agr.sp$AGS[germany_agr.sp$AGS %in% berid] <- "11000"

germany_agr.df <- unique(data.frame(germany_agr.sp[, c("AGS")]) )
germany_agr.sp <- unionSpatialPolygons(germany_agr.sp, germany_agr.sp$AGS)

row.names(germany_agr.df) <- germany_agr.df$AGS
germany_agr.sp <- SpatialPolygonsDataFrame(germany_agr.sp, germany_agr.df)


### Overall borders (https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-1-000-000-ebenen-stand-01-01-vg1000-ebenen-01-01.html)
ger.sp <- readOGR(dsn = "./vg1000_2019-01-01.utm32s.shape.ebenen/vg1000_ebenen",
                  layer = "VG1000_STA") 
proj4string(ger.sp)
ger.sp <- gBuffer(ger.sp[1, ], byid = FALSE, width = 500)




#################################
### Read census 2011 age data ###
#################################
# https://www.regionalstatistik.de/genesis/online/data;sid=FE602C43F3C36B8FA84AAA849E84D6E7.reg2?operation=abruftabelleAbrufen&selectionname=12111-04-01-4&levelindex=0&levelid=1585048532267&index=14



zensus_age.df <- read.csv("12111-04-01-4.csv", skip = 9, header = FALSE, sep = ";",
                          quote = "\"", stringsAsFactors = FALSE, nrows = 54876, na.strings = "-",
                          colClasses = c("character"))
names(zensus_age.df) <- c("AGS", "Landkreis", "key", "value", "value_male", "value_female")


### Clean key
zensus_age.df$key <- gsub(" bis unter ", "_", zensus_age.df$key)
zensus_age.df$key <- gsub(" Jahre| und mehr|unter ", "", zensus_age.df$key)
zensus_age.df$key <- gsub(" Jahr", "", zensus_age.df$key)
zensus_age.df$key <- paste0("age_", zensus_age.df$key)


### Spread key
zensus_age.df <- tidyr::spread(zensus_age.df[, c("AGS", "Landkreis", "key", "value")], 
                               key = "key", value = "value", convert = TRUE)

### Numeric
oo <- grep("age_", names(zensus_age.df))
zensus_age.df[, oo] <- apply(zensus_age.df[, oo], 2, FUN = as.numeric)


### Combine to RKI age groups
age_gr<- list(0:5, 5:15, 15:35, 35:60, 60:80, 80:100)

for(i in age_gr){
  i <- unlist(i)
  names <- paste0("age_", i[1:(length(i)-1)], "_", i[2:(length(i))] )
  if(i[1] == 0){
    names[1] <- "age_1"
  }
  if(i[length(i)] == 100){
    names <- c(names, "age_100")
  }
  newnames <- paste0("age_gr_", min(i), "_", (max(i) - 1))
  
  zensus_age.df[, newnames] <- rowSums(zensus_age.df[, names], na.rm = TRUE)

}


### Reduce information
oo <- grep("age_gr_|age_Insgesamt", names(zensus_age.df))
zensus_age.df <- zensus_age.df[, c(1:2, oo)]





#######################################
### Aggregate age to covid counties ###
#######################################

### Change AGS
zensus_age.df$AGS[zensus_age.df$AGS == "02"] <- "02000"
zensus_age.df$AGS[zensus_age.df$AGS == "11"] <- "11000"

# Goettingen doubled
zensus_age.df <- zensus_age.df[-which(zensus_age.df$AGS == "03159"), ]
zensus_age.df$AGS[zensus_age.df$AGS == "03152"] <- "03159"

# # Berlin
# oo <- which(nchar(zensus_age.df$AGS) == 8 & substr(zensus_age.df$AGS, 1, 2) == "11")
# zensus_age.df$AGS[oo] <- substr(zensus_age.df$AGS[oo], 1, 5)


### Aggregate over Change in Mecklenburg Vorpommern
p <- as.character(c(13071:13076))
zensus_age.df <- zensus_age.df[-which(zensus_age.df$AGS %in% p), ]

p <- c("13056", "13055", "13002", "13052")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Mecklenburgische Seenplatte"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13071"

p <- c("13053", "13051")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Rostock"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13072"

p <- c("13005", "13057", "13061")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Vorpommern-Rügen"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13073"

p <- c("13006", "13058")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Nordwestmecklenburg"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13074"

p <- c("13001", "13059", "13062")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Vorpommern-Greifswald"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13075"

p <- c("13054", "13060")
zensus_age.df$Landkreis[zensus_age.df$AGS %in% p] <- "LK Ludwigslust-Parchim"
zensus_age.df$AGS[zensus_age.df$AGS %in% p] <- "13076"

# Aggregate
oo <- which(names(zensus_age.df) %in% c("AGS", "Landkreis"))
zensus_age.df <- aggregate(zensus_age.df[, -oo],
                           by = zensus_age.df[, oo],
                           FUN = sum)




###################
### Proportions ###
###################

vars <- names(zensus_age.df)[grep("age_gr", names(zensus_age.df))]

for(i in vars){
  newname <- paste0("prop_", i)
  zensus_age.df[, newname] <- zensus_age.df[, i] / zensus_age.df$age_Insgesamt
}


##################
### Merge data ###
##################

# ### Test
# test <- germany_covid_agr.df[germany_covid_agr.df$date == max(germany_covid_agr.df$date),]
# test2 <- merge(test, zensus_age.df, 
#                 by = "AGS", all.x = TRUE, all.y = TRUE)
# View(test2[which(is.na(test2$date)), c("AGS", "Landkreis.y", "age_Insgesamt")])


### Merge
germany_covid_cens.df <- merge(germany_covid_agr.df, zensus_age.df, 
                               by = "AGS", all.x = TRUE, all.y = FALSE)






##################################
### Plot population age groups ###
##################################

### Most recent date

germany_covid_recent.df <- germany_covid_cens.df[germany_covid_cens.df$date == max(germany_covid_cens.df$date),]


### Cutoff points
mycut <- function(x, n = 15, t = 5, p = 0.9, start = NULL){
  x <- x[which(!is.na(x))]
  
  p1 <- seq(p, 0.99, length.out = t)
  suppressWarnings(c1 <- quantileCuts(x, params = p1))
  
  if(is.null(start)){
    p_st <- (1:((n - t) - 1))/(n - t)
    start <- p_st[1]
  }
  
  p2 <- seq(start, p - (p1[2] - p1[1]) * 1.5, length.out = (n - t))
  suppressWarnings(c2 <- quantileCuts(x, params = p2))
  
  
  return(c(c2, c1))
}




### Variables and names Names


agevars <- c("prop_age_gr_0_4" ,"prop_age_gr_5_14", "prop_age_gr_15_34", "prop_age_gr_35_59", "prop_age_gr_60_79", "prop_age_gr_80_99")
names(agevars) <- c("% 0-4 years old", "% 5-14 years old", "% 15-34 years old",
                    "% 35-59 years old", "% 60-79 years old", "% 80 years and older")
summary(rowSums(germany_covid_recent.df[, agevars]))

agevars_abs <- c("age_gr_0_4" ,"age_gr_5_14", "age_gr_15_34", "age_gr_35_59", "age_gr_60_79", "age_gr_80_99")
names(agevars_abs) <- c("0-4 years old", "5-14 years old", "15-34 years old",
                    "35-59 years old", "60-79 years old", "80 years and older")
sum(rowSums(germany_covid_recent.df[, agevars_abs]))
sum(germany_covid_recent.df$age_Insgesamt)



### Plot share ###




# Variables and combined vector
allvalues <- unlist(germany_covid_recent.df[, agevars])

# Colours
cols2 <- inferno(16, begin = 0.1, end = 1, direction = -1)
cols2 <- c(cols2)


# All cases identical cut-off
vacant.shades1 <- auto.shading(allvalues[which(allvalues > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols2, digits = 4)
cuts <- mycut(allvalues[which(allvalues > 0)], 
              n = 15, t = 3, p = 0.92)
vacant.shades1$breaks <- c(cuts)

vacant.shades1


### Merge shape and data
tmp.spdf <- merge(germany_agr.sp, germany_covid_recent.df, by = "AGS")

### Overall percentage

t1 <- weighted.mean(data.frame(tmp.spdf)[, agevars[1]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t2 <- weighted.mean(data.frame(tmp.spdf)[, agevars[2]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t3 <- weighted.mean(data.frame(tmp.spdf)[, agevars[3]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t4 <- weighted.mean(data.frame(tmp.spdf)[, agevars[4]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t5 <- weighted.mean(data.frame(tmp.spdf)[, agevars[5]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t6 <- weighted.mean(data.frame(tmp.spdf)[, agevars[6]], w = tmp.spdf$age_Insgesamt, na.rm = TRUE)

t <- round(c(t1, t2, t3, t4, t5, t6), 2)



### Plot

png(file = paste0("../03_Output/", "Germany_age_", "2011", ".png"), width = 14, height = 10, 
    units = "in", bg = "white", family = "CM Roman", res = 400)
par(mar = c(2, 0, 3, 6))
par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))


# Loop over age groups

for(i in 1:length(agevars)){
  choropleth(tmp.spdf, data.frame(tmp.spdf)[, agevars[i]], shading = vacant.shades1, border = NA,
             main = paste0(names(agevars)[i], ": ", t[i]), cex.main = 1.6)
  
  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.10), y2, cex = 1.2, vacant.shades1, title = "Percent",
               border = NA, fmt = "%.2f", x.intersp = 1.2)
  par(xpd = FALSE)
}



### Outer label
mtext(paste0("Distribution of age groups percent (2011)"), outer = TRUE, cex = 1.5, line = 1)

mtext("Data source: Zensus 2011, https://www.regionalstatistik.de/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = -2)
mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = 0)

dev.off()





### Plot absolute values ###


# Variables and combined vector
allvalues <- unlist(germany_covid_recent.df[, agevars_abs]) / 1000

# Colours
cols2 <- inferno(16, begin = 0.1, end = 1, direction = -1)
cols2 <- c(cols2)


# All cases identical cut-off
vacant.shades1 <- auto.shading(allvalues[which(allvalues > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols2, digits = 4)
cuts <- mycut(allvalues[which(allvalues > 0)], 
              n = 15, t = 3, p = 0.92)
vacant.shades1$breaks <- c(cuts)


### Merge shape and data
tmp.spdf <- merge(germany_agr.sp, germany_covid_recent.df, by = "AGS")

### Overall percentage

t1 <- sum(data.frame(tmp.spdf)[, agevars_abs[1]] / 1000, na.rm = TRUE)
t2 <- sum(data.frame(tmp.spdf)[, agevars_abs[2]] / 1000, na.rm = TRUE)
t3 <- sum(data.frame(tmp.spdf)[, agevars_abs[3]] / 1000, na.rm = TRUE)
t4 <- sum(data.frame(tmp.spdf)[, agevars_abs[4]] / 1000, na.rm = TRUE)
t5 <- sum(data.frame(tmp.spdf)[, agevars_abs[5]] / 1000, na.rm = TRUE)
t6 <- sum(data.frame(tmp.spdf)[, agevars_abs[6]] / 1000, na.rm = TRUE)

t <- round(c(t1, t2, t3, t4, t5, t6), 0)



### Plot

png(file = paste0("../03_Output/", "Germany_age_abs", "2011", ".png"), width = 14, height = 10, 
    units = "in", bg = "white", family = "CM Roman", res = 400)
par(mar = c(2, 0, 3, 6))
par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))


# Loop over age groups

for(i in 1:length(agevars_abs)){
  choropleth(tmp.spdf, data.frame(tmp.spdf)[, agevars_abs[i]] / 1000, 
             shading = vacant.shades1, border = NA,
             main = paste0(names(agevars_abs)[i], ": ", t[i]), cex.main = 1.6)
  
  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.13), y2, cex = 1.2, vacant.shades1, title = "Absolute (in 1,000) ",
               border = NA, fmt = "%.0f")
  par(xpd = FALSE)
}



### Outer label
mtext(paste0("Absolute number in each age group (2011)"), outer = TRUE, cex = 1.5, line = 1)

mtext("Data source: Zensus 2011, https://www.regionalstatistik.de/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = -2)
mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = 0)

dev.off()








#############################################################
### Difference in shares from each group (60- 79 and 80+) ###
#############################################################


### Compute proportion of positive tests
covvars <- c("sum_cases_00_04", "sum_cases_05_14", "sum_cases_15_34",
            "sum_cases_35_59", "sum_cases_60_79", "sum_cases_80")

propcovvars <- NULL
for(i in covvars){
  newname <- paste0("prop_", gsub("sum_", "", i))
  germany_covid_recent.df[, newname] <- germany_covid_recent.df[, i] / germany_covid_recent.df$sum_cases
  germany_covid_recent.df[, newname][is.infinite(germany_covid_recent.df[, newname])] <- NA
  propcovvars <- c(propcovvars, newname)
}

summary(germany_covid_recent.df[, propcovvars], na.rm = TRUE)


names(propcovvars) <- c("% Cases  0-4 years old", "% Cases 5-14 years old", "% Cases 15-34 years old", 
                 "% Cases 35-59 years old", "% Cases 60-79 years old",  "% Cases 80 years and older") 
names(agevars) <- gsub("%", "% Population", names(agevars), fixed = TRUE)


### Difference between shares
germany_covid_recent.df$diff_0_4 <-   germany_covid_recent.df$prop_cases_00_04 - germany_covid_recent.df$prop_age_gr_0_4   
germany_covid_recent.df$diff_5_14 <-  germany_covid_recent.df$prop_cases_05_14 - germany_covid_recent.df$prop_age_gr_5_14  
germany_covid_recent.df$diff_15_34 <- germany_covid_recent.df$prop_cases_15_34 - germany_covid_recent.df$prop_age_gr_15_34 
germany_covid_recent.df$diff_35_59 <- germany_covid_recent.df$prop_cases_35_59 - germany_covid_recent.df$prop_age_gr_35_59 
germany_covid_recent.df$diff_60_79 <- germany_covid_recent.df$prop_cases_60_79 - germany_covid_recent.df$prop_age_gr_60_79 
germany_covid_recent.df$diff_80 <-    germany_covid_recent.df$prop_cases_80 - germany_covid_recent.df$prop_age_gr_80_99 

summary(germany_covid_recent.df$diff_60_79)
summary(germany_covid_recent.df$diff_80)

diffvars <- c("diff_0_4", "diff_5_14", "diff_15_34", "diff_35_59", "diff_60_79", "diff_80")
names(diffvars) <- c("Difference 0-4 years old", "Difference 5-14 years old", 
                     "Difference 15-34 years old", "Difference 35-59 years old", 
                     "Difference 60-79 years old", "Difference 80 years and older")



### Summary by BL

for(i in diffvars){
  cat(paste("\n\n\n\n", i, "\n\n"))
  print(by(germany_covid_recent.df[, i], germany_covid_recent.df$BL, 
           FUN = function(x) summary(x)))
}




############################################
### Plot age groups group 60- 79 and 80+ ###
############################################


# Variables and combined vector
allvalues1 <- unlist(germany_covid_recent.df[, c(agevars[5:6])]) * 100
allvalues2 <- unlist(germany_covid_recent.df[, c(propcovvars[5:6])]) * 100
allvalues3 <- unlist(germany_covid_recent.df[, diffvars[5:6]]) * 100


# Colours
cols1 <- viridis(16, begin = 0, end = 1, direction = -1)
cols1 <- c(cols1)
cols2 <- inferno(16, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#D3D3D3", cols2)


cols3 <- c(viridis(11, begin = 0.3, end = 0.9, direction = 1),
           brewer.pal(5, "OrRd"))




# Percent of positive and population
vacant.shades1 <- auto.shading(allvalues1[which(allvalues1 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols1, digits = 4)
cuts <- mycut(allvalues1[which(allvalues1 > 0)], 
              n = 15, t = 3, p = 0.92)
vacant.shades1$breaks <- c(cuts)

vacant.shades1


# Percent of positive and population
vacant.shades2 <- auto.shading(allvalues2[which(allvalues2 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols2, digits = 4)
cuts2 <- mycut(allvalues2[which(allvalues2 > 0)], 
              n = 15, t = 2, p = 0.90)
vacant.shades2$breaks <- c(0.0000001, cuts2)

vacant.shades2



# Difference
vacant.shades3 <- auto.shading(allvalues3[which(!is.na(allvalues3))], 
                               cutter = quantileCuts, 
                               n = 16, col = cols3, digits = 4)
cuts3 <- mycut(allvalues3, 
              n = 15, t = 3, p = 0.92)
cuts3_rev <- rev(-mycut(-allvalues3, 
                   n = 15, t = 3, p = 0.92))
cuts3_mid <- quantileCuts(allvalues3[which(!is.na(allvalues3))] , 
                          params = c(seq(0.14, 0.86, length.out = 9)))
vacant.shades3$breaks <- c(cuts3_rev[1:3], cuts3_mid, cuts3[13:15])

vacant.shades3



### Merge shape and data
tmp.spdf <- merge(germany_agr.sp, germany_covid_recent.df, by = "AGS")

### Overall percentage

t1p <- weighted.mean(data.frame(tmp.spdf)[, agevars[5]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t2p <- weighted.mean(data.frame(tmp.spdf)[, agevars[6]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)

t1c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[5]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)
t2c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[6]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)

t1d <- t1c - t1p
t2d <- t2c - t2p

t <- round(c(t1p, t1c, t1d, t2p,  t2c, t2d), 2)

leg <- c("Per 100 inhabitants", "Per 100 cases", "Percentage-points",
         "Per 100 inhabitants", "Per 100 cases", "Percentage-points")
under <- c("under", "", "under",
           "under", "", "under")

shades <- list(vacant.shades1, vacant.shades2, vacant.shades3,
               vacant.shades1, vacant.shades2, vacant.shades3)


### Plot
date <- unique(tmp.spdf$date)

png(file = paste0("../03_Output/", "Germany_age60_differences_", date, ".png"), width = 14, height = 10, 
    units = "in", bg = "white", family = "CM Roman", res = 400)
par(mar = c(2, 0, 3, 6))
par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))


# Loop over age groups

plotvars <- c(agevars[5], propcovvars[5], diffvars[5],
              agevars[6], propcovvars[6], diffvars[6])
for(i in 1:length(plotvars)){
  choropleth(tmp.spdf, data.frame(tmp.spdf)[, plotvars[i]] * 100, shading = shades[[i]], border = NA,
             main = paste0(names(plotvars)[i], ": ", t[i]), cex.main = 1.6)
  
  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.13), y2, cex = 1.2, shades[[i]], title = leg[i],
               border = NA, fmt = "%.2f", x.intersp = 1.5, under = under[i])
  par(xpd = FALSE)
}



### Outer label
mtext(paste0("Relative COVID-19 cases among age groups 60-79 and 80+: ", date), outer = TRUE, cex = 1.5, line = 1)

mtext("Data source: RKI, https://npgeo-corona-npgeo-de.hub.arcgis.com/ & Zensus 2011, https://www.regionalstatistik.de/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = -2)
mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = 0)

dev.off()







#####################################
### Plot age groups group 15 - 49 ###
#####################################


# Variables and combined vector
allvalues1 <- unlist(germany_covid_recent.df[, c(agevars[3:4])]) * 100
allvalues2 <- unlist(germany_covid_recent.df[, c(propcovvars[3:4])]) * 100
allvalues3 <- unlist(germany_covid_recent.df[, diffvars[3:4]]) * 100


# Colours
cols1 <- viridis(16, begin = 0, end = 1, direction = -1)
cols1 <- c(cols1)
cols2 <- inferno(16, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#D3D3D3", cols2)


cols3 <- c(viridis(6, begin = 0.3, end = 0.9, direction = 1),
           brewer.pal(9, "OrRd"), "#361223")




# Percent of positive and population
vacant.shades1 <- auto.shading(allvalues1[which(allvalues1 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols1, digits = 4)
cuts <- mycut(allvalues1[which(allvalues1 > 0)], 
              n = 15, t = 3, p = 0.92)
vacant.shades1$breaks <- c(cuts)

vacant.shades1


# Percent of positive and population
vacant.shades2 <- auto.shading(allvalues2[which(allvalues2 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols2, digits = 4)
cuts2 <- mycut(allvalues2[which(allvalues2 > 0)], 
               n = 15, t = 2, p = 0.90)
vacant.shades2$breaks <- c(0.0000001, cuts2)

vacant.shades2



# Difference
vacant.shades3 <- auto.shading(allvalues3[which(!is.na(allvalues3))], 
                               cutter = quantileCuts, 
                               n = 16, col = cols3, digits = 4)
cuts3 <- mycut(allvalues3, 
               n = 15, t = 3, p = 0.92)
cuts3_rev <- rev(-mycut(-allvalues3, 
                        n = 15, t = 3, p = 0.92))
cuts3_mid <- quantileCuts(allvalues3[which(!is.na(allvalues3))] , 
                          params = c(seq(0.14, 0.86, length.out = 9)))
vacant.shades3$breaks <- c(cuts3_rev[1:3], cuts3_mid, cuts3[13:15])

vacant.shades3



### Merge shape and data
tmp.spdf <- merge(germany_agr.sp, germany_covid_recent.df, by = "AGS")

### Overall percentage

t1p <- weighted.mean(data.frame(tmp.spdf)[, agevars[3]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t2p <- weighted.mean(data.frame(tmp.spdf)[, agevars[4]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)

t1c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[3]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)
t2c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[4]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)

t1d <- t1c - t1p
t2d <- t2c - t2p

t <- round(c(t1p, t1c, t1d, t2p,  t2c, t2d), 2)

leg <- c("Per 100 inhabitants", "Per 100 cases", "Percentage-points",
         "Per 100 inhabitants", "Per 100 cases", "Percentage-points")
under <- c("under", "", "under",
           "under", "", "under")

shades <- list(vacant.shades1, vacant.shades2, vacant.shades3,
               vacant.shades1, vacant.shades2, vacant.shades3)


### Plot
date <- unique(tmp.spdf$date)

png(file = paste0("../03_Output/", "Germany_age15_differences_", date, ".png"), width = 14, height = 10, 
    units = "in", bg = "white", family = "CM Roman", res = 400)
par(mar = c(2, 0, 3, 6))
par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))


# Loop over age groups

plotvars <- c(agevars[3], propcovvars[3], diffvars[3],
              agevars[4], propcovvars[4], diffvars[4])
for(i in 1:length(plotvars)){
  choropleth(tmp.spdf, data.frame(tmp.spdf)[, plotvars[i]] * 100, shading = shades[[i]], border = NA,
             main = paste0(names(plotvars)[i], ": ", t[i]), cex.main = 1.6)
  
  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.13), y2, cex = 1.2, shades[[i]], title = leg[i],
               border = NA, fmt = "%.2f", x.intersp = 1.5, under = under[i])
  par(xpd = FALSE)
}



### Outer label
mtext(paste0("Relative COVID-19 cases among age groups 15-34 and 35-59: ", date), outer = TRUE, cex = 1.5, line = 1)

mtext("Data source: RKI, https://npgeo-corona-npgeo-de.hub.arcgis.com/ & Zensus 2011, https://www.regionalstatistik.de/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = -2)
mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = 0)

dev.off()







######################################
### Plot age groups group below 15 ###
######################################


# Variables and combined vector
allvalues1 <- unlist(germany_covid_recent.df[, c(agevars[1:2])]) * 100
allvalues2 <- unlist(germany_covid_recent.df[, c(propcovvars[1:2])]) * 100
allvalues3 <- unlist(germany_covid_recent.df[, diffvars[1:2]]) * 100


# Colours
cols1 <- viridis(16, begin = 0, end = 1, direction = -1)
cols1 <- c(cols1)
cols2 <- inferno(16, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#D3D3D3", cols2)


cols3 <- c(viridis(14, begin = 0.3, end = 0.9, direction = 1),
           brewer.pal(3, "OrRd")[-2])




# Percent of positive and population
vacant.shades1 <- auto.shading(allvalues1[which(allvalues1 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols1, digits = 4)
cuts <- mycut(allvalues1[which(allvalues1 > 0)], 
              n = 15, t = 3, p = 0.92)
vacant.shades1$breaks <- c(cuts)

vacant.shades1


# Percent of positive and population
vacant.shades2 <- auto.shading(allvalues2[which(allvalues2 > 0)], 
                               cutter = quantileCuts,  
                               n = 16, col = cols2, digits = 4)
cuts2 <- mycut(allvalues2[which(allvalues2 > 0)], 
               n = 15, t = 2, p = 0.90)
vacant.shades2$breaks <- c(0.0000001, cuts2)

vacant.shades2



# Difference
vacant.shades3 <- auto.shading(allvalues3[which(!is.na(allvalues3))], 
                               cutter = quantileCuts, 
                               n = 16, col = cols3, digits = 4)
cuts3 <- mycut(allvalues3, 
               n = 15, t = 3, p = 0.92)
cuts3_rev <- rev(-mycut(-allvalues3, 
                        n = 15, t = 3, p = 0.92))
cuts3_mid <- quantileCuts(allvalues3[which(!is.na(allvalues3))] , 
                          params = c(seq(0.14, 0.86, length.out = 9)))
vacant.shades3$breaks <- c(cuts3_rev[1:3], cuts3_mid, cuts3[13:15])

vacant.shades3



### Merge shape and data
tmp.spdf <- merge(germany_agr.sp, germany_covid_recent.df, by = "AGS")

### Overall percentage

t1p <- weighted.mean(data.frame(tmp.spdf)[, agevars[1]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)
t2p <- weighted.mean(data.frame(tmp.spdf)[, agevars[2]] * 100, w = tmp.spdf$age_Insgesamt, na.rm = TRUE)

t1c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[1]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)
t2c <- weighted.mean(data.frame(tmp.spdf)[, propcovvars[2]] * 100, w = tmp.spdf$sum_cases, na.rm = TRUE)

t1d <- t1c - t1p
t2d <- t2c - t2p

t <- round(c(t1p, t1c, t1d, t2p,  t2c, t2d), 2)

leg <- c("Per 100 inhabitants", "Per 100 cases", "Percentage-points",
         "Per 100 inhabitants", "Per 100 cases", "Percentage-points")
under <- c("under", "", "under",
           "under", "", "under")

shades <- list(vacant.shades1, vacant.shades2, vacant.shades3,
               vacant.shades1, vacant.shades2, vacant.shades3)


### Plot
date <- unique(tmp.spdf$date)

png(file = paste0("../03_Output/", "Germany_age0_differences_", date, ".png"), width = 14, height = 10, 
    units = "in", bg = "white", family = "CM Roman", res = 400)
par(mar = c(2, 0, 3, 6))
par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))


# Loop over age groups

plotvars <- c(agevars[1], propcovvars[1], diffvars[1],
              agevars[2], propcovvars[2], diffvars[2])
for(i in 1:length(plotvars)){
  choropleth(tmp.spdf, data.frame(tmp.spdf)[, plotvars[i]] * 100, shading = shades[[i]], border = NA,
             main = paste0(names(plotvars)[i], ": ", t[i]), cex.main = 1.6)
  
  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.13), y2, cex = 1.2, shades[[i]], title = leg[i],
               border = NA, fmt = "%.2f", x.intersp = 1.5, under = under[i])
  par(xpd = FALSE)
}



### Outer label
mtext(paste0("Relative COVID-19 cases among age groups 0-4 and 5-15: ", date), outer = TRUE, cex = 1.5, line = 1)

mtext("Data source: RKI, https://npgeo-corona-npgeo-de.hub.arcgis.com/ & Zensus 2011, https://www.regionalstatistik.de/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = -2)
mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
      cex = 1.2, side = 1, adj = 1, line = 0)

dev.off()



