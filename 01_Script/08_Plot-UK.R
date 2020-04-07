
#### Plot RKI Covid Germany ####
#### Tobias Ruettenauer ####
#### 2020/ 03 / 21 ####

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




#################
### Load data ###
#################
### UK data prepared by emmadoughty (https://github.com/emmadoughty/Daily_COVID-19) 


ukd <- "https://raw.githubusercontent.com/emmadoughty/Daily_COVID-19/master/Data/cases_by_area.csv"


### UK data
uk.df <- read.table(ukd,
                    header = TRUE, sep = ",", na.strings = "NA",
                    quote = "\"")


####################
### Prepare data ###
####################

### Date
names(uk.df)[1] <- "date"
uk.df$date <- as.Date(uk.df$date, format = "%d/%m/%Y")

### Keep only UTLA in England
uk.df <- uk.df[uk.df$type == "UTLA" & uk.df$country == "England", ]

### Keep only cases since 07/03
uk.df <- uk.df[uk.df$date >= "2020-03-07", ]

### rename id
uk.df$ctyua19cd <- uk.df$GSS_CD

# Cases at recent date
sum(uk.df$confirm[uk.df$date == max(uk.df$date)])



### New cases for eacht day
uk.df <- uk.df[order(uk.df$area, uk.df$date), ]

uk.df$daily_cases <- ave(uk.df$confirm, 
                       by = uk.df$GSS_CD,
                       FUN = function(x) x - dplyr::lag(x))


######################
### Load shapefile ###
######################
# using https://geoportal.statistics.gov.uk/datasets/counties-and-unitary-authorities-december-2019-boundaries-uk-buc/data?page=6

uk.sp <- readOGR(dsn = "./Counties_and_Unitary_Authorities_December_2019_Boundaries_UK_BUC",
                        layer = "Counties_and_Unitary_Authorities_December_2019_Boundaries_UK_BUC") 
proj4string(uk.sp)

### Keep only England
uk.sp$country <- substr(as.character(uk.sp$ctyua19cd), 1, 1)
uk.sp <- uk.sp[uk.sp$country == "E", ]


### Drop Isles of Scilly (E06000053), as part of Cornwell in covid data
uk.sp <- uk.sp[-which(uk.sp$ctyua19cd == "E06000053"),]


### Overall borders (save and load, as this is really slow)
uk_border.sp <- unionSpatialPolygons(uk.sp, ID = rep(1, nrow(uk.sp)))
uk_border.sp <- gBuffer(uk_border.sp, byid = FALSE, width = 500)




### Aggregate to data units




###-----------------------------------###
###               Plot                ###
###-----------------------------------###



#######################
### listwise object ###
#######################

### Nb
uk.nb <- poly2nb(uk.sp)

# Replace empty nb with neares neigbour
uk.1nb <- knn2nb(knearneigh(coordinates(uk.sp), k = 1))
uk.nb[which(card(uk.nb) == 0)] <- uk.1nb[which(card(uk.nb) == 0)]

### Listwise object
uk.lw <- nb2listw(uk.nb, style = "W")




##############################
### Define cuts and colors ###
##############################

### Colours
cols1 <- viridis(15, direction = -1)
cols1 <- c("#FAF0E6", cols1)

cols2 <- inferno(14, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#B0C4DE", "#FAF0E6", cols2 )



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
vacant.shades1 <- auto.shading(uk.df$confirm[which(uk.df$confirm > 0)], 
                               cutter = quantileCuts,  
                               n = 15, col = cols1, digits = 4)
cuts <- mycut(uk.df$confirm[which(uk.df$confirm > 0)], 
              n = 14, t = 4, p = 0.92)
vacant.shades1$breaks <- c(0.0000000000001, cuts)

# New cases
vacant.shades2 <- auto.shading(uk.df$daily_cases[which(uk.df$daily_cases > 0)], 
                               cutter = quantileCuts,  
                               n = 12, col = cols2, digits = 4)
cuts2 <- mycut(uk.df$daily_cases[which(uk.df$daily_cases > 0)], 
              n = 13, t = 3, p = 0.92)
vacant.shades2$breaks <- c(0, 0.0000000000001, cuts2)




#########################
### Plot cases by day ###
#########################

dates <- unique(uk.df$date)


for(i in dates){

  ### Merge shape and data
  tmp.df <- uk.df[uk.df$date == i, ]

  tmp.spdf <- merge(uk.sp, tmp.df,
                    by = "ctyua19cd")

  tc <- sum(tmp.df$confirm, na.rm = TRUE)
  nc <- sum(tmp.df$daily_cases, na.rm = TRUE)

  j <- as.character(as.Date(i, origin = "1970-01-01"))
  
  ### Fill zeros (city of London)
  tmp.spdf$confirm[is.na(tmp.spdf$confirm)] <- 0
  tmp.spdf$daily_cases[is.na(tmp.spdf$daily_cases)] <- 0
  

  ### Autocorrelation
  if(i != dates[1]){
    # moran_all <- moran.test(tmp.spdf$confirm, uk.lw)
    moran_new <- moran.test(tmp.spdf$daily_cases, uk.lw)
  }
  moran_all <- moran.test(tmp.spdf$confirm, uk.lw)
  # moran_new <- moran.test(tmp.spdf$daily_cases, uk.lw)



  ### Plot

  png(file = paste0("../03_Output/", "England_cases_", j, ".png"), width = 12, height = 7,
      units = "in", bg = "white", family = "CM Roman", res = 600)
  par(mar=c(2, 0, 3, 6))
  par(mfrow=c(1, 2), oma = c(2, 0, 1.5, 0))


  #### Total cases

  choropleth(tmp.spdf, tmp.spdf$confirm, shading = vacant.shades1, border = NA,
             main = paste0("Total cases: ", tc), cex.main = 1.5)

  plot(tmp.spdf, border = ggplot2::alpha("grey70", 0.5), lwd = 0.5, add = T)
  plot(uk_border.sp, border = "orange1", lwd = 1, add = T)


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
  plot(uk_border.sp, border = "orange1", lwd = 1, add = T)


  # Coordinates of window
  x1 <- par()$usr[1]
  x2 <- par()$usr[2]
  y1 <- par()$usr[3]
  y2 <- par()$usr[4]
  r <- x2 - x1

  # Legend
  par(xpd = NA)
  choro.legend((x2 - r*0.10), y2, cex = 1, vacant.shades2, title = "New cases ",
               border = NA, fmt = "%.0f")
  par(xpd = FALSE)

  if(i == dates[1]){
    mtext(paste0("Moran's I:     "), outer = FALSE, cex = 1.5, side = 1)
  }else{
    mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  }
  # mtext(paste0("Moran's I: ", round(moran_new$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)

  ### Outer label
  mtext(paste0("COVID-19: ", j), outer = TRUE, cex = 1.5)

  mtext("Data source: Emma Doughty, https://github.com/emmadoughty/Daily_COVID-19", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 0)
  mtext("Code and processed data: https://github.com/ruettenauer/COVID-19-maps/", outer = TRUE,
        cex = 0.8, side = 1, adj = 1, line = 1)

  dev.off()


}



######################
### Animate as gif ###
######################

# All pictures
files <- paste0("../03_Output/", "England_cases_", dates, ".png")

# Repeat last picture
files <- c(files, files[length(files)], files[length(files)], files[length(files)], files[length(files)]) 

# Import files
img <- lapply(files, FUN = function(x) image_read(x))

# Resize
img2 <- lapply(img, FUN = function(x) image_resize(x, "1260x735"))

# Animate
img2 <- image_join(img2)

gif <- image_animate(img2, fps = 2)

image_write(gif, "../03_Output/England_covid19.gif")





