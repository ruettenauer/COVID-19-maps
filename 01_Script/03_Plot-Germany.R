
#### Read Covid Data from John Hopkins ####
#### Tobias Ruettenauer ####
#### 2020/ 01 / 31 ####

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

### WHO Data directory (https://github.com/CSSEGISandData/COVID-19)
whod <- "C:/work/Forschung/Daten/COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

### Italy data directory (https://github.com/pcm-dpc/COVID-19)
itad <- "C:/work/Forschung/Daten/COVID-19-ita/dati-province/"


#################
### Load data ###
#################


load("Germany_covid.RData")


######################
### Load shapefile ###
######################
# using RKI shapefile (https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/917fc37a709542548cc3be077a786c17_0)

germany.sp <- readOGR(dsn = "./RKI_Corona_Landkreise",
                        layer = "RKI_Corona_Landkreise") 
proj4string(germany.sp)
germany.sp <- spTransform(germany.sp, CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

### Replace AGS with RS for Berlin
germany.sp$AGS <- as.character(germany.sp$AGS)
germany.sp$RS <- as.character(germany.sp$RS)
germany.sp$AGS[which(is.na(germany.sp$AGS))] <- germany.sp$RS[which(is.na(germany.sp$AGS))]

### Overall borders (https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-1-000-000-ebenen-stand-01-01-vg1000-ebenen-01-01.html)
ger.sp <- readOGR(dsn = "./vg1000_2019-01-01.utm32s.shape.ebenen/vg1000_ebenen",
                      layer = "VG1000_STA") 
proj4string(ger.sp)
ger.sp <- gBuffer(ger.sp[1, ], byid = FALSE, width = 500)




##########################
### Gen balanced panel ###
##########################

### Drop NA Landkreis
germany_long.df <- germany_long.df[which(germany_long.df$AGS != "0-1"),]

### Get all counties and dates since 2020-02-25
germany_lk.df <- data.frame(germany.sp)

mind <- as.Date("2020-02-25", format = "%Y-%m-%d")
maxd <- max(germany_long.df$date)
dates <- seq.Date(mind, maxd, 1)

### Gen full df
germany_covid.df <- data.frame(AGS = rep(germany_lk.df$AGS, times = length(dates)),
                              date = rep(dates, each = nrow(germany_lk.df)))



### Add columns
vars <- c("AGS", "SDV_RS", "GEN", "NUTS", "EWZ", "KFL", "BL", "BL_ID")
germany_covid.df <- merge(germany_covid.df, germany_lk.df[, vars],
                         by = "AGS", all.x = TRUE)

### Merge daily reports
drop <- c("id_bl", "Bundesland")
germany_covid.df <- merge(germany_covid.df, germany_long.df[, -c(which(names(germany_long.df) %in% drop))], 
                         by = c("AGS", "date"), all.x = TRUE, all.y = TRUE)

### Fill cumulative sums

#Repeat last function (credits to https://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value)
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}

# get cumsums
cumvars <- which(grepl("sum_", names(germany_covid.df)))

# fill last NA obs
germany_covid.df <- germany_covid.df[order(germany_covid.df$AGS, germany_covid.df$date), ]
germany_covid.df[, cumvars] <- apply(germany_covid.df[, cumvars], 2, 
                                 FUN = function(x) ave(x, 
                                                       by = germany_covid.df$AGS,
                                                       FUN = function(z) repeat_last(z)))


### Fill zeros for report variables
repvars <- which(!names(germany_covid.df) %in% c(names(germany_lk.df), "date", "Landkreis", "Meldedatum"))
germany_covid.df[, repvars][is.na(germany_covid.df[, repvars])] <- 0


# Total cases 2020-03-20 (as of 2020-03-21)
sum(germany_covid.df$sum_cases[germany_covid.df$date == "2020-03-20"], na.rm = TRUE)


### Save
save(germany_covid.df, file = "Germany_covid_full.RData")


###-----------------------------------###
###               Plot                ###
###-----------------------------------###



#######################
### listwise object ###
#######################

de.nb <- poly2nb(germany.sp)
de.lw <- nb2listw(de.nb, style = "W")




##############################
### Define cuts and colors ###
##############################

cols1 <- viridis(15, direction = -1)
cols1 <- c("#FAF0E6", cols1)

cols2 <- inferno(15, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#FAF0E6", cols2 )

# Cutoff for two quantile distributions
midpoint <- quantile(germany_covid.df$sum_cases[which(germany_covid.df$sum_cases > 0)], 
                     na.rm = TRUE, probs = c(0.9, 0.95, 0.99))

midpoint2 <- quantile(germany_covid.df$daily_cases[which(germany_covid.df$sum_cases > 0)], 
                     na.rm = TRUE, probs = c(0.9, 0.95, 0.99))


# Total cases
vacant.shades1 <- auto.shading(germany_covid.df$sum_cases[which(germany_covid.df$sum_cases > 0)], n = 16, cutter = quantileCuts, 
                               col = cols1, digits = 4)

cuts1 <- quantileCuts(germany_covid.df$sum_cases[which(germany_covid.df$sum_cases > 0 &
                                                         germany_covid.df$sum_cases <= midpoint[2])], 9)

cuts2 <- quantileCuts(germany_covid.df$sum_cases[which(germany_covid.df$sum_cases > midpoint[1])], 4)


vacant.shades1$breaks <- c(0.0000000000001, cuts1, midpoint[1], cuts2, max(germany_covid.df$sum_cases)/3)

# cuts <- scales::trans_breaks("log2", function(x) 2 ^ x, 16)(germany_covid.df$sum_cases)
# 
# cutsvacant.shades1$breaks <- c(0.0000000000001, cuts)




# New cases
vacant.shades2 <- auto.shading(germany_covid.df$daily_cases[which(germany_covid.df$daily_cases > 0)], n = 15, 
                               cutter = quantileCuts, 
                               col = cols2, digits = 4)

if(length(vacant.shades2$breaks) < 15){
  vacant.shades2$breaks <- c(vacant.shades2$breaks, max(germany_covid.df$daily_cases)/2)
  
  cols2 <- inferno((length(vacant.shades2$breaks) + 2), begin = 0.1, end = 1, direction = -1)
  cols2 <- c("#FAF0E6", cols2 )
  vacant.shades2$cols <- cols2
}


vacant.shades2$breaks <- c(0.0000000000001, vacant.shades2$breaks)



############################################
### Plot cases and negative cases by day ###
############################################

dates <- unique(germany_covid.df$date)

### Most recent date seems incomplete -> omit
dates <- dates[which(dates >= "2020-02-25" & dates < max(dates))]

for(i in dates){
  
  ### Merge shape and data
  tmp.df <- germany_covid.df[germany_covid.df$date == i, ]
  
  tmp.spdf <- merge(germany.sp, tmp.df,
                    by = "AGS")
  
  tc <- sum(tmp.df$sum_cases, na.rm = TRUE)
  nc <- sum(tmp.df$daily_cases, na.rm = TRUE)
  
  j <- as.character(as.Date(i, origin = "1970-01-01"))
  
  ### Autocorrelation
  # if(i != dates[1]){
  #   moran_all <- moran.test(tmp.spdf$sum_cases, de.lw)
  #   moran_new <- moran.test(tmp.spdf$daily_cases, de.lw)
  # }
  moran_all <- moran.test(tmp.spdf$sum_cases, de.lw)
  moran_new <- moran.test(tmp.spdf$daily_cases, de.lw)

  
  
  ### Plot
  
  png(file = paste0("../03_Output/", "Germany_cases_", j, ".png"), width = 12, height = 7, 
      units = "in", bg = "white", family = "CM Roman", res = 600)
  par(mar=c(2, 0, 3, 6))
  par(mfrow=c(1, 2), oma = c(2, 0, 1.5, 0))
  
  
  #### Total cases
  
  choropleth(tmp.spdf, tmp.spdf$sum_cases, shading = vacant.shades1, border = NA,
             main = paste0("Total cases: ", tc), cex.main = 1.5)
  
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
  choro.legend((x2 - r*0.10), y2, cex = 1, vacant.shades1, title = "Total cases ",
               border = NA, fmt = "%.0f", under = "")
  par(xpd = FALSE)
  
  # if(i == dates[1]){
  #   mtext(paste0("Moran's I:     "), outer = FALSE, cex = 1.5, side = 1)
  # }else{
  #   mtext(paste0("Moran's I: ", round(moran_all$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  # }
  mtext(paste0("Moran's I: ", round(moran_all$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  #### New cases
  
  choropleth(tmp.spdf, tmp.spdf$daily_cases, shading = vacant.shades2, border = NA,
             main = paste0("New cases: ", nc), cex.main = 1.5)
  
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
  choro.legend((x2 - r*0.10), y2, cex = 1, vacant.shades2, title = "New cases ",
               border = NA, fmt = "%.0f", under = "")
  par(xpd = FALSE)
  
  # if(i == dates[1]){
  #   mtext(paste0("Moran's I:     "), outer = FALSE, cex = 1.5, side = 1)
  # }else{
  #   mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  # }
  mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  ### Outer label
  mtext(paste0("COVID-19: ", j), outer = TRUE, cex = 1.5)
  
  mtext("Data source: RKI, https://npgeo-corona-npgeo-de.hub.arcgis.com/", outer = TRUE, 
        cex = 0.8, side = 1, adj = 1, line = 0)
  mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE, 
        cex = 0.8, side = 1, adj = 1, line = 1)
  
  dev.off()
  
  
}



######################
### Animate as gif ###
######################


# files <- list.files(path = "../03_Output/", pattern = "Germany_cases_*", full.names = T)
files <- paste0("../03_Output/", "Germany_cases_", dates, ".png")


img <- lapply(files, FUN = function(x) image_read(x))

# img2 <- lapply(img, FUN = function(x) image_resize(x, "840x490"))

img2 <- lapply(img, FUN = function(x) image_resize(x, "1260x735"))


img2 <- image_join(img2)

gif <- image_animate(img2, fps = 1)

image_write(gif, "../03_Output/Germany_covid19.gif")










