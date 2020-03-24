
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


###################################
### Read and prepare INKAR data ###
###################################
# Source: https://www.inkar.de/ (manually downloaded)

inkardir <- "./INKAR/"

i <- dir(inkardir)

j <- i[grepl("INKAR", i)]
j <- paste0(inkardir, j)

c <- 1

for(k in j){
  
  # Clean header
  header <- as.vector(t(read.table(k, nrows = 1, sep = ";")[1, ]))
  header <- gsub(" ", "", header)

  # remove remaining non ascii characters
  header <- iconv(header, "latin1", "ASCII", sub="")
  header <- make.names(header)
  header[1] <- "kennziffer"
  
  # # Combine with second row header (year)
  # header2 <- as.vector(t(read.table(k, skip=1, nrows = 1, sep = ";")[1,]))
  # header3 <- paste(header, header2, sep="_")
  # header3 <- gsub("_NA", "", header3)
  
  # Input and rename data
  data <- read.csv(k, skip = 2, header = FALSE, sep = ";",
                   quote = "\"", dec = ",", stringsAsFactors = FALSE,
                   colClasses = c(V1 = "character"))
  names(data) <- header
  
  # Correct character vars (containing thousands separator)
  vars <- which(sapply(data, is.character))
  vars <- vars[-which(vars %in% c(1:3))]
  for(l in vars){
    data[, l] <- gsub("\\.", "", data[, l])
    data[, l] <- gsub("\\,", ".", data[, l])
    data[, l] <- as.numeric(data[, l])
  }
  
  # Combine 
  if(k != j[1]){
    inkar_2017.df <- merge(inkar_2017.df, data, all.x = TRUE,
                           by = c("kennziffer", "Raumeinheit", "Aggregat"))
  }else{
    inkar_2017.df <- data
  }
  

  c <- c + 1
  
}



### Add meta data
# https://www.inkar.de/documents/Referenz%20Gemeinden-GVB-Kreise_NUTS.xlsx

meta.df <- read_xlsx("Referenz Gemeinden-GVB-Kreise_NUTS.xlsx", sheet = "KRS")

# Correct names
names <- names(meta.df)
oo <- grep("..", names, fixed = TRUE)
names[oo] <- paste0(names[oo - 1], "name")
names(meta.df) <- names

# Drop explanation line
meta.df <- meta.df[-1, ]

# Are and inhabitants numeric

meta.df[, c("fl17", "bev17")] <- apply(meta.df[, c("fl17", "bev17")], 2, 
                                       FUN = function(x) as.numeric(x))

# Make identical id
ids <- as.numeric(meta.df$krs17)
newids <- format(ids, scientific = FALSE)
newids <- gsub(" ", "", newids)
newids[which(ids <  10000000)] <- paste0("0", newids[which(ids <  10000000)], sep = "")
newids <- substr(newids, 1, 5)
meta.df$kennziffer <- newids

# Merge
inkar_2017.df <- merge(meta.df, inkar_2017.df, by = "kennziffer", all.x = TRUE, all.y = TRUE)



##################
### Merge data ###
##################

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


### merge with INKAR

germany_covid_inkar.df <- merge(germany_covid_agr.df, inkar_2017.df,
                                by.x = "AGS", by.y = "kennziffer",
                                all.x = TRUE, all.y = TRUE)


### Save combined data
save(germany_covid_inkar.df, file = "Germany_covid_inkar.RData")






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




