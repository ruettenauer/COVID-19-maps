
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
library(ggplot2)
library(scales)

library(extrafont)
loadfonts()


### Working Directory
setwd("C:/work/Forschung/Covid-19/02_Data")




#################
### Load data ###
#################

dit <- "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv"

### Italian total data
italy.df <- read.table(dit,
                       header = TRUE, sep = ",", na.strings = "",
                       quote = "\"")
# Format date
italy.df$date <- as.Date(substr(italy.df$data, 1, 10), "%Y-%m-%d")



####################
### Prepare data ###
####################

italy.df <- italy.df[order(italy.df$date), ]

### Daily tests
italy.df$new_tests <- italy.df$tamponi - dplyr::lag(italy.df$tamponi)

### fraction positive
italy.df$new_positive_per <- italy.df$nuovi_positivi / italy.df$new_tests * 100
italy.df <- italy.df[!is.na(italy.df$new_positive_per), ]



########################
### Barchart per day ###
########################

max1 <- max(italy.df$new_positive_per)
max2 <- max(italy.df$new_tests)
r <- max2/max1

zp1 <- ggplot(italy.df, aes(x = date, y = new_positive_per))
zp1 <- zp1 + geom_bar(stat = "identity", aes(y = new_positive_per, fill = "Fraction positive tests"))
zp1 <- zp1 + geom_text(aes(label = paste0(sprintf("%.1f", round(new_positive_per, 1)), "")), 
                       vjust = -0.5,  size = 4.5, colour = "steelblue", family = "CM Roman",
                       check_overlap = TRUE)
zp1 <- zp1 + geom_line(aes(x = date, y = new_tests/r), colour = alpha("deepskyblue", 0.8), lwd = 1.6)
zp1 <- zp1 + geom_line(aes(x = date, y = new_tests/r), colour = alpha("deepskyblue1", 0.5), lwd = 1.3)
zp1 <- zp1 + geom_line(aes(x = date, y = new_tests/r), colour = alpha("deepskyblue2", 0.5), lwd = 1.0)
zp1 <- zp1 + geom_line(aes(x = date, y = new_tests/r), colour = alpha("deepskyblue3", 0.5), lwd = 0.5)
#zp1 <- zp1 + geom_line(aes(x = date, y = new_tests/r), colour = "deepskyblue4", lwd = 1)
zp1 <- zp1 + geom_point(aes(y = new_tests/r, shape = "Absolute number of tests"), 
                        size = 3, colour = "orangered2")
zp1 <- zp1 + scale_shape_manual(values = 18)
zp1 <- zp1 + scale_color_manual(values = "orangered2")
zp1 <- zp1 + scale_fill_manual(values = "grey35")
zp1 <- zp1 + scale_x_date(breaks = pretty_breaks(16), labels = date_format("%d.%m"))
zp1 <- zp1 + scale_y_continuous(sec.axis = sec_axis(~.*r, name = "Absolute number of tests per day"))
zp1 <- zp1 + xlab("Time") + ylab("% positive tests per day") 
zp1 <- zp1 + ggtitle("COVID-19 Italy: Fraction of positive tests") 
zp1 <- zp1 + theme_bw()
zp1 <- zp1 + theme(text = element_text(family = "CM Roman", size = 20),
                   axis.text.y = element_text(colour = "black"),
                   axis.text.x = element_text(colour = "black", size = 12),
                   axis.title.x = element_text(colour = "black",  margin = margin(t = 10, r = 0 , b = 0, l = 0)),
                   axis.title.y = element_text(colour = "black",  margin = margin(t = 0, r = 15 , b = 0, l = 0)),
                   axis.title.y.right = element_text(colour = "black",  margin = margin(t = 0, r = 0 , b = 0, l = 15)),
                   plot.title = element_text(hjust = 0.5),
                   plot.margin = unit(c(1,1,1,1), "lines"),
                   plot.caption = element_text(size = 12, margin = margin(t = 15), lineheight = 0.5),
                   legend.key = element_blank(), legend.title = element_blank(),
                   legend.text = element_text(size = 12),
                   legend.position = c(0.01,0.95), legend.justification = 'left',
                   legend.background = element_blank(),
                   legend.box.background = element_rect(colour = "black"),
                   legend.spacing.y = unit(-0.1, "cm"))
zp1 <- zp1 + labs(caption = "Data source: Dipartimento della Protezione Civile, https://github.com/pcm-dpc/COVID-19 \n 
                     Code and processed data: https://github.com/ruettenauer/COVID-19-maps/")

zp1


# Export
png(file = paste0("../03_Output/", "Italy_positive-fraction", ".png"), width = 16, height = 9,
    units = "in", bg = "white", family = "CM Roman", res = 600)
par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
zp1
dev.off()


