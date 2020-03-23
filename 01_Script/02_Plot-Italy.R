
#### Plot italy covid ####
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




#################
### Load data ###
#################


load("Italy_covid.RData")


######################
### Load shapefile ###
######################

provinces.sp <- readOGR(dsn = "./Limiti01012019_g/ProvCM01012019_g",
                        layer = "ProvCM01012019_g_WGS84") 

provinces.sp <- spTransform(provinces.sp, CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


# Overall borders
italy.sp <- readOGR(dsn = "./Limiti01012019_g/RipGeo01012019_g",
                        layer = "RipGeo01012019_g_WGS84") 

italy.sp <- spTransform(italy.sp, CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


italy.sp <- gBuffer(italy.sp, byid = FALSE, width = 500)



###-----------------------------------###
###               Plot                ###
###-----------------------------------###



#######################
### listwise object ###
#######################

it.nb <- poly2nb(provinces.sp)
it.lw <- nb2listw(it.nb, style = "W")




############################################
### Define cuts and colors ###
############################################

cols1 <- viridis(15, direction = -1)
cols1 <- c("#FAF0E6", cols1)

cols2 <- inferno(14, begin = 0.1, end = 1, direction = -1)
cols2 <- c("#B0C4DE", "#FAF0E6", cols2 )

# Total cases
vacant.shades1 <- auto.shading(italy.df$totale_casi[!is.na(italy.df$totale_casi)], n = 16, cutter = quantileCuts, 
                               col = cols1, digits = 4)

med <- median(italy.df$totale_casi)
cuts1 <- quantileCuts(italy.df$totale_casi[!is.na(italy.df$totale_casi)
                                           & italy.df$totale_casi > 0
                                           & italy.df$totale_casi <= 200], 8)

cuts2 <- quantileCuts(italy.df$totale_casi[!is.na(italy.df$totale_casi)
                                           & italy.df$totale_casi > 0
                                           & italy.df$totale_casi > 200], 7)

vacant.shades1$breaks <- c(0.0000000000001, cuts1, cuts2)

# cuts <- scales::trans_breaks("log2", function(x) 2 ^ x, 16)(italy.df$totale_casi)
# 
# cutsvacant.shades1$breaks <- c(0.0000000000001, cuts)




# New cases
vacant.shades2 <- auto.shading(italy.df$new_cases[which(italy.df$new_cases > 0)], n = 15, cutter = quantileCuts, 
                               col = cols2, digits = 4)

vacant.shades2$breaks <- c(0, 0.0000000000001, vacant.shades2$breaks)



############################################
### Plot cases and negative cases by day ###
############################################

dates <- unique(italy.df$date)

for(i in dates){
  
  ### Merge shape and data
  tmp.df <- italy.df[italy.df$date == i,]
  
  tmp.spdf <- merge(provinces.sp, tmp.df,
                    by = "COD_PROV")
  
  tc <- sum(tmp.df$totale_casi, na.rm = TRUE)
  nc <- sum(tmp.df$new_cases, na.rm = TRUE)
  
  j <- as.character(as.Date(i, origin = "1970-01-01"))
  
  ### Autocorrelation
  if(i != dates[1]){
    moran_all <- moran.test(tmp.spdf$totale_casi, it.lw)
    moran_new <- moran.test(tmp.spdf$new_cases, it.lw)
  }

  
  
  ### Plot
  
  png(file = paste0("../03_Output/", "Italy_cases_", j, ".png"), width = 12, height = 7, 
      units = "in", bg = "white", family = "CM Roman", res = 600)
  par(mar=c(2, 0, 3, 6))
  par(mfrow=c(1, 2), oma = c(2, 0, 1.5, 0))
  
  
  #### Total cases
  
  choropleth(tmp.spdf, tmp.spdf$totale_casi, shading = vacant.shades1, border = NA,
             main = paste0("Total cases: ", tc), cex.main = 1.5)
  
  plot(tmp.spdf, border = "grey70", lwd = 1, add = T)
  plot(italy.sp, border = "orange1", lwd = 1, add = T)
  
  
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
  
  if(i == dates[1]){
    mtext(paste0("Moran's I:     "), outer = FALSE, cex = 1.5, side = 1)
  }else{
    mtext(paste0("Moran's I: ", round(moran_all$estimate[1], 2)), outer = FALSE, cex = 1.5, side = 1)
  }
  
  
  #### New cases
  
  choropleth(tmp.spdf, tmp.spdf$new_cases, shading = vacant.shades2, border = NA,
             main = paste0("New cases: ", nc), cex.main = 1.5)
  
  plot(tmp.spdf, border = "grey70", lwd = 1, add = T)
  plot(italy.sp, border = "orange1", lwd = 1, add = T)
  
  
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
  
  
  ### Outer label
  mtext(paste0("COVID-19: ", j), outer = TRUE, cex = 1.5)
  
  mtext("Data source: https://github.com/pcm-dpc/COVID-19", outer = TRUE, 
        cex = 1, side = 1, adj = 1)
  
  dev.off()
  
  
}



######################
### Animate as gif ###
######################


files <- list.files(path = "../03_Output/", pattern = "Italy_cases_*", full.names = T)

img <- lapply(files, FUN = function(x) image_read(x))

# img2 <- lapply(img, FUN = function(x) image_resize(x, "840x490"))

img2 <- lapply(img, FUN = function(x) image_resize(x, "1260x735"))


img2 <- image_join(img2)

gif <- image_animate(img2, fps = 1)

image_write(gif, "../03_Output/Italy_covid19.gif")










