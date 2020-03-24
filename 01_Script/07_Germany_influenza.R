
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



### Weekly covid data
week_covid.df <- read.csv("survstat_covid_2020-03-23/Data.csv", skip = 1, header = TRUE, sep = "\t",
                          quote = "\"", stringsAsFactors = FALSE, skipNul = TRUE)

names(week_covid.df)[1] <- "county"
names(week_covid.df) <- gsub("X2020.", "covid_20.", names(week_covid.df))


### Weekly influenza data
week_influenza.df <- read.csv("survstat_influenza_2020-03-23/Data.csv", skip = 1, header = TRUE, sep = "\t",
                          quote = "\"", stringsAsFactors = FALSE, skipNul = TRUE)

names(week_influenza.df)[1] <- "county"
names(week_influenza.df) <- gsub("X2017.|X2018.", "infl_1718.", names(week_influenza.df))





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


### Overall borders (https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-1-000-000-ebenen-stand-01-01-vg1000-ebenen-01-01.html)
ger.sp <- readOGR(dsn = "./vg1000_2019-01-01.utm32s.shape.ebenen/vg1000_ebenen",
                  layer = "VG1000_STA") 
proj4string(ger.sp)
ger.sp <- gBuffer(ger.sp[1, ], byid = FALSE, width = 500)



############################
### Combine IDs and fill ###
############################

meta.df <- data.frame(germany.sp)


### Clean German special characters
meta.df$county <- gsub("[.]{2,}", ".", make.names(iconv(meta.df$county, "latin1", "ASCII", sub="")))
week_influenza.df$county <- gsub("[.]{2,}", ".", make.names(iconv(week_influenza.df$county, "latin1", "ASCII", sub="")))
week_covid.df$county <- gsub("[.]{2,}", ".", make.names(iconv(week_covid.df$county, "latin1", "ASCII", sub="")))


# Clean manually
meta.df$county[which(meta.df$county == "LK.Neustadt.a.d.Aisch.Bad.Windsheim")] <- "LK.Neustadt.Aisch.Bad.Windsheim"
meta.df$county[which(meta.df$county == "StadtRegion.Aachen")] <- "StdteRegion.Aachen"


### Merge with ids
week_influenza.df <- merge(meta.df, week_influenza.df, by = "county", all.x = TRUE)
week_covid.df <- merge(meta.df, week_covid.df, by = "county", all.x = TRUE)


### fill NAs
oo <- grep(".KW", names(week_influenza.df))
week_influenza.df[, oo][is.na(week_influenza.df[, oo])] <- 0

oo <- grep(".KW", names(week_covid.df))
week_covid.df[, oo][is.na(week_covid.df[, oo])] <- 0




####################################
### Reshape to long format files ###
####################################

oo <- grep(".KW", names(week_influenza.df))
week_influenza.df <- reshape(week_influenza.df,  idvar = "OBJECTID",
                             varying = oo, sep = ".",
                             direction = "long")

oo <- grep(".KW", names(week_covid.df))
week_covid.df <- reshape(week_covid.df, idvar = "OBJECTID",
                             varying = oo, sep = ".",
                             direction = "long")


### Combine

cases_weekly.df <- merge(week_influenza.df, week_covid.df[, c("OBJECTID", "time", "covid_20")],
                         by = c("OBJECTID", "time"), all.x = TRUE, all.y = TRUE)



### Gen season week
w <- as.character(c(27:52, 1:26))
w[which(nchar(w) < 2)] <- paste0("0", w[which(nchar(w) < 2)])
sw <- paste0("KW", w)
tmp <- data.frame(time = sw, week = c(1: length(sw)))

cases_weekly.df <- merge(cases_weekly.df, tmp, by = "time",
                         all.x = TRUE)

sort <- c("OBJECTID", "time", "week")
cases_weekly.df <- cases_weekly.df[, c(sort, names(cases_weekly.df)[-which(names(cases_weekly.df) %in% sort)])]



### Fill non-existing weeks if smaller / equal than max sw
m1 <- max(c(1:52)[which(sw %in% week_influenza.df$time)])
m2 <- max(c(1:52)[which(sw %in% week_covid.df$time)])


cases_weekly.df$infl_1718[which(is.na(cases_weekly.df$infl_1718) & cases_weekly.df$week <= m1)] <- 0
cases_weekly.df$covid_20[which(is.na(cases_weekly.df$covid_20) & cases_weekly.df$week <= m2)] <- 0


### Generate cumulative numbers
cases_weekly.df <- cases_weekly.df[order(cases_weekly.df$OBJECTID, cases_weekly.df$week), ] 
  
cases_weekly.df$sum_infl_1718 <- ave(cases_weekly.df$infl_1718,
                                     by = cases_weekly.df$OBJECTID,
                                     FUN = function(x) cumsum(x))

cases_weekly.df$sum_covid_20<- ave(cases_weekly.df$covid_20,
                                     by = cases_weekly.df$OBJECTID,
                                     FUN = function(x) cumsum(x))




###-----------------------------------###
###               Plot                ###
###-----------------------------------###



#######################
### listwise object ###
#######################

de.nb <- poly2nb(germany.sp)
de.lw <- nb2listw(de.nb, style = "W")



##############################
### Plot influenza by week ###
##############################

by(cases_weekly.df$sum_infl_1718, cases_weekly.df$week, sum)

# Plot weeks 20 - 43

data.df <- cases_weekly.df[which(cases_weekly.df$week >= 20 & cases_weekly.df$week <= 43),]


##############################
### Define cuts and colors ###
##############################

### Colours
cols1 <- viridis(15, direction = -1)
cols1 <- c("#FAF0E6", cols1)

cols2 <- inferno(12, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#FAF0E6", cols2 )



### Cutoff points
mycut <- function(x, n = 15, t = 5, p = 0.9, start = NULL){
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


# Total cases
vacant.shades1 <- auto.shading(data.df$sum_infl_1718[which(data.df$sum_infl_1718 > 0)], 
                               cutter = quantileCuts,  
                               n = 15, col = cols1, digits = 4)
cuts <- mycut(data.df$sum_infl_1718[which(data.df$sum_infl_1718 > 0)], 
              n = 15, t = 3, p = 0.92, start = 0.05)
vacant.shades1$breaks <- c(0.0000000000001, cuts[-1])

# New cases
vacant.shades2 <- auto.shading(data.df$infl_1718[which(data.df$infl_1718 > 0)], 
                               cutter = quantileCuts,  
                               n = 12, col = cols2, digits = 4)
cuts2 <- mycut(data.df$infl_1718[which(data.df$infl_1718 > 0)], 
               n = 11, t = 3, p = 0.94, start = 0.15)
vacant.shades2$breaks <- c(0.0000000000001, cuts2[-1], max(data.df$infl_1718)/2)





#########################
### Plot cases by day ###
#########################


for(i in unique(data.df$week)){
  
  ### Merge shape and data
  tmp.df <- data.df[data.df$week == i, ]
  
  tmp.spdf <- merge(germany.sp, tmp.df,
                    by = "RS")
  
  tc <- sum(tmp.df$sum_infl_1718, na.rm = TRUE)
  nc <- sum(tmp.df$infl_1718, na.rm = TRUE)
  
  j <- 2017
  if(i > 	26){j <- j + 1}
  date <- paste0(i, " (", gsub("KW", "CW", unique(tmp.df$time)), " ", j, ")" )
  
  ### Autocorrelation
  # if(i != dates[1]){
  #   moran_all <- moran.test(tmp.spdf$sum_cases, de.lw)
  #   moran_new <- moran.test(tmp.spdf$daily_cases, de.lw)
  # }
  moran_all <- moran.test(tmp.spdf$sum_infl_1718, de.lw)
  moran_new <- moran.test(tmp.spdf$infl_1718, de.lw)
  
  
  
  ### Plot
  
  png(file = paste0("../03_Output/", "Germany_influenza17_w", i, ".png"), width = 12, height = 7,
      units = "in", bg = "white", family = "CM Roman", res = 400)
  par(mar=c(2, 0, 3, 6))
  par(mfrow=c(1, 2), oma = c(2, 0, 1.5, 0))
  
  
  #### Total cases
  
  choropleth(tmp.spdf, tmp.spdf$sum_infl_1718, shading = vacant.shades1, border = NA,
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
  
  # Morans I
  mtext(paste0("Moran's I: ", round(moran_all$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  
  #### New cases
  choropleth(tmp.spdf, tmp.spdf$infl_1718, shading = vacant.shades2, border = NA,
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
  
  # Morans I
  mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  ### Outer label
  mtext(paste0("Influenza 17/18: season week ", date), outer = TRUE, cex = 1.5)
  
  mtext("Data source: Robert Koch-Institut, https://survstat.rki.de", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 0)
  mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 1)
  
  dev.off()
  
  
}








###-----------------------------------###
###  Plot Compare Influenza covid     ###
###-----------------------------------###


### Since case 100 ###

sum <- ave(germany_covid.df$sum_cases,
            by = germany_covid.df$date,
            FUN = function(x) sum(x))

data_covid.df <- germany_covid.df[which(sum > 100),]

sum <- ave(cases_weekly.df$sum_infl_1718,
           by = cases_weekly.df$week,
           FUN = function(x) sum(x))

data_infl.df <- cases_weekly.df[which(sum > 100),]


### Cut to same time length

tp <- length(unique(data_covid.df$date))

data_infl.df$tnr <- ave(data_infl.df$sum_infl_1718,
           by = data_infl.df$OBJECTID,
           FUN = function(x) c(1:length(x)))
data_infl.df <- data_infl.df[which(data_infl.df$tnr <= tp), ]


##############################
### Define cuts and colors ###
##############################

### Colours
cols1 <- viridis(15, direction = -1)
cols1 <- c("#FAF0E6", cols1)

cols2 <- inferno(15, begin = 0, end = 1, direction = -1)
cols2 <- c("#FAF0E6", cols2 )



### Cutoff points
mycut <- function(x, n = 15, t = 5, p = 0.9, start = NULL){
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


# Total cases influenza
vacant.shades1 <- auto.shading(data_infl.df$sum_infl_1718[which(data_infl.df$sum_infl_1718 > 0)], 
                               cutter = quantileCuts,  
                               n = 15, col = cols1, digits = 4)
cuts <- mycut(data_infl.df$sum_infl_1718[which(data_infl.df$sum_infl_1718 > 0)], 
              n = 15, t = 3, p = 0.94, start = 0.15)
vacant.shades1$breaks <- c(0.0000000000001, cuts[-1])


# Total cases covid
vacant.shades2 <- auto.shading(data_covid.df$sum_cases[which(data_covid.df$sum_cases > 0)], 
                               cutter = quantileCuts,  
                               n = 15, col = cols2, digits = 4)
cuts2 <- mycut(data_covid.df$sum_cases[which(data_covid.df$sum_cases > 0)], 
              n = 14, t = 3, p = 0.92, start = 0.10)
vacant.shades2$breaks <- c(0.0000000000001, cuts2)






################################
### Plot cases by week / day ###
################################

c <- c(1: length(unique(data_covid.df$date)))
weeks <- unique(data_infl.df$week)
days <- unique(data_covid.df$date)


for(i in c){
  
  ### Merge shape and data
  tmp_covid.df <- data_covid.df[data_covid.df$date == days[i], ]
  tmp_infl.df <- data_infl.df[data_infl.df$week == weeks[i], ]
  
  
  tmp_covid.spdf <- merge(germany.sp, tmp_covid.df,
                          by = "AGS")
  tmp_infl.spdf <- merge(germany.sp, tmp_infl.df,
                          by = "AGS")
  
  ti <- sum(tmp_infl.spdf$sum_infl_1718, na.rm = TRUE)
  tc <- sum(tmp_covid.spdf$sum_cases, na.rm = TRUE)
  

  ### Autocorrelation
  # if(i != dates[1]){
  #   moran_all <- moran.test(tmp.spdf$sum_cases, de.lw)
  #   moran_new <- moran.test(tmp.spdf$daily_cases, de.lw)
  # }
  moran_all <- moran.test(tmp_infl.spdf$sum_infl_1718, de.lw)
  moran_new <- moran.test(tmp_covid.spdf$sum_cases, de.lw)
  
  
  
  ### Plot
  
  png(file = paste0("../03_Output/", "Germany_influenza-covid_", i, ".png"), width = 12, height = 7,
      units = "in", bg = "white", family = "CM Roman", res = 400)
  par(mar=c(2, 0, 4, 6))
  par(mfrow=c(1, 2), oma = c(2, 0, 1.5, 0))
  
  
  #### Total cases
  
  choropleth(tmp_infl.spdf, tmp_infl.spdf$sum_infl_1718, shading = vacant.shades1, border = NA,
             main = paste0("Influenza cases: ", ti, "\n", "Week: ", i ), cex.main = 1.5)
  
  plot(tmp_infl.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
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
  
  # Morans I
  mtext(paste0("Moran's I: ", round(moran_all$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  
  #### New cases
  choropleth(tmp_covid.spdf, tmp_covid.spdf$sum_cases, shading = vacant.shades2, border = NA,
             main = paste0("COVID-19 cases: ", tc, "\n", "Day: ", i ), cex.main = 1.5)
  
  plot(tmp_covid.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(ger.sp, border = "orange1", lwd = 1, add = T)
  
  
  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1
  
  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.10), y2, cex = 1, vacant.shades2, title = "Total cases ",
               border = NA, fmt = "%.0f", under = "")
  par(xpd = FALSE)
  
  # Morans I
  mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  
  ### Outer label
  mtext(paste0("Influenza 17/18 in weeks and COVID-19 in days since case 100"), outer = TRUE, cex = 1.5)
  
  mtext("Data source: Robert Koch-Institut, https://survstat.rki.de & https://npgeo-corona-npgeo-de.hub.arcgis.com/", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 0)
  mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 1)
  
  dev.off()
  
  
}





######################
### Animate as gif ###
######################


# Omit last two day (because of reporting lag)
f <- c[1:(length(c) - 2)]

# Repeat last picture
f <- c(f, f[length(f)], f[length(f)]) 

# Import files
files <- paste0("../03_Output/", "Germany_influenza-covid_", f, ".png")

img <- lapply(files, FUN = function(x) image_read(x))

# Resize
img2 <- lapply(img, FUN = function(x) image_resize(x, "1260x735"))

# Animate
img2 <- image_join(img2)

gif <- image_animate(img2, fps = 1)

image_write(gif, "../03_Output/Germany_influenza_covid.gif")




