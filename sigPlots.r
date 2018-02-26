## Load libraries

library(tidyverse)
library(devtools)
library(RSQLite)
library(grid)
library(gridExtra)

library(rworldmap)
library(ggmap)
library(plotly)
library(htmlwidgets)

Sys.setenv(TZ = 'GMT')

########
### Import datasets into list
########

# Define Project ID
projectID <- 109

source("loadMotusData.r")
tagList <- loadMotusData(projectID)
bandSite_only <- TRUE

# Left join receivers with tagList to ensure any NAs are corrected if possible 
# Only overwrites NAs or 0s in original detection file
if (bandSite_only) {
  fixed_tagList <- tagList %>%
    filter(site %in% c(bpsites, sealsites)) %>%
    mutate(lat = ifelse(is.na(lat), latitude, ifelse(lat==0,latitude, lat)), 
           lon = ifelse(is.na(lon), longitude, ifelse(lon==0,longitude, lon)))
} else {
  fixed_tagList <- tagList %>%
    filter(!site %in% shitList) %>%
    left_join(receivers, by = 'recv') %>%
    mutate(lat = ifelse(is.na(lat), latitude, ifelse(lat==0,latitude, lat)), 
           lon = ifelse(is.na(lon), longitude, ifelse(lon==0,longitude, lon)), 
           site = ifelse(site %in% bpsites, 'Bon Portage', ifelse(site %in% sealsites, 'Seal Island', ifelse(is.na(site), deploymentName, site))))
}


flights.raw <- read_csv("Project109-sigplots-Flights.csv") 
flights.raw <- flights.raw %>%
  mutate(depRunID = ifelse(is.na(depRunID), max(flights.raw[flights.raw$markerNumber == markerNumber, 5:16]), depRunID))

flight.meta <- flights.raw %>% select(bandsite, markerNumber, code, depRunID)

flights.df <- fixed_tagList[fixed_tagList$runID %in% as.numeric(unlist(flights.raw[5:16])),]

flights <- flights.df %>%
  group_by(runID) %>%
  summarise(ts = ts[which.max(sig)], hitID = hitID[which.max(sig)],sig = max(sig), markerNumber = markerNumber[1], age = age[1])

depFlights <- flights %>%
  left_join(flight.meta, by = 'markerNumber') %>%
  group_by(markerNumber) %>%
  summarise(depRunID = ifelse(is.na(depRunID), runID[which.max(ts)], max(depRunID))[1], 
            ts = ts[runID==depRunID], hitID = hitID[runID==depRunID], sig = sig[runID==depRunID], runID = runID[runID==depRunID], 
            bandsite = bandsite[1], age = age[1], code = code[1]) %>%
  select(-depRunID)


####################################
## Create a list of plots of signal strength/time/ant for each tag at and each site 
## At the end there is the option to save the plots for later

# Make a grid of plots for each tower visted by each tag
for (i in unique(fixed_tagList$markerNumber)) {
  print(i)
#  i <- '1861-32538'
  ## Select a single tag and convert to tibble
  markerSigs <- fixed_tagList %>% filter(markerNumber == i) %>% mutate(site = ifelse(is.na(site), ifelse(is.na(recv), 'UNK', recv), site)) %>%
    mutate(ts = as.POSIXct(ts, origin = '1970-01-01'), antBearing = as.factor(ifelse(is.na(antBearing), 'UNK', antBearing)))
#  markerSigs <- filter(markerSigs, day(ts) == 17)
  markerFlights <- flights %>% filter(markerNumber == i)
  markerDepFlight <- depFlights %>% filter(markerNumber == i)
  
  ## Check if it there is any data
  if (nrow(markerSigs) > 0 & !is.na(i)) {
    
  sunRiseSetTimes <- tibble(ts = days <- seq(min(markerSigs$ts)-days(1), max(markerSigs$ts), 24*3600),
                            lat = rep(43.4674, length(ts)),
                            lon = rep(-65.7489, length(ts)))
  
  sunRiseSetTimes <- sunRiseSetTimes %>% sunRiseSet(lat = "lat", lon = "lon", ts = "ts") %>%
    mutate(date = as.Date(ts)) %>%
    select(-lat, -lon, -ts)
  
  attr(sunRiseSetTimes$sunset, 'tzone') <- "America/Halifax"
  attr(sunRiseSetTimes$sunrise, 'tzone') <- "America/Halifax"
  
  bandingDates <- rawBandingData$DATE.BANDED[rawBandingData$markerNumber == i]+hours(4)  # Add 4 hours to set to GMT
  
#  bandingDates <- sort(bandingDates, decreasing = F)[-1]
  
  if (nrow(markerSigs %>% filter(year(ts) == year(min(bandingDates)))) > 0) {
    bandingDates <- bandingDates[year(bandingDates)==min(year(bandingDates))]
    markerSigs <- markerSigs %>% filter(year(ts) == year(min(bandingDates)))
  } else {
    bandingDates <- bandingDates[year(bandingDates)==max(year(bandingDates))]
    markerSigs <- markerSigs %>% filter(year(ts) == year(max(bandingDates)))
  }
  
  latestDetection <- ifelse(max(bandingDates) > max(markerSigs$ts), max(bandingDates), max(markerSigs$ts))
#  earliestDetection <- min(bandingDates)
  earliestDetection <- min(bandingDates)

  ## Create a list of plots, one for each site
  markerSigs %>% 
      ggplot() + # Setting the date to something easy to interpret!
      scale_x_datetime(date_labels = "%d/%m/%y %H:%M:%S", limits = c(earliestDetection, latestDetection)) + # Format the x-axis so it displays the full date
      xlab("Timestamp (D/M/Y H:M:S)") +
      ggtitle(paste0(i,' - ', unique(markerSigs$age))) +
      geom_vline(xintercept = min(bandingDates), color = 'black', size = 2)+
      geom_vline(xintercept = bandingDates, color = 'black', size = 1)+
      geom_vline(xintercept = sunRiseSetTimes$sunrise, color = 'red')+
      geom_vline(xintercept = sunRiseSetTimes$sunset, color = 'blue')+
      geom_point(aes(as.POSIXct(ts, origin = '1970-01-01'), sig, color = antBearing)) + 
      geom_point(data = markerFlights, aes(ts, sig), size = 4, shape = 1, stroke = 2) + 
      geom_point(data = markerDepFlight, aes(ts, sig), size = 3, shape = 4, stroke = 2) + 
      guides(color = FALSE) +
      facet_grid(site+antBearing~.)
      ggsave(paste0('proj#',projectID,'-tagPlot#', i ,'.png'), scale = 5)
   
    p <- markerSigs %>% 
      group_by(runID) %>%
      summarise(sig = max(sig), ts = median(ts), antBearing = antBearing[1], site = site[1], age = age[1]) %>%
      mutate(antBearing = as.factor(ifelse(is.na(antBearing), 'UNK', antBearing)), ts = as.POSIXct(ts), site = ifelse(is.na(site), 'UNK', site)) %>% # Change NA's into "UNK"
      ggplot() + # Setting the date to something easy to interpret!
    #  scale_x_datetime(date_labels = "%d/%m/%y %H:%M:%S", limits = c(min(bandingDates), latestDetection)) + # Format the x-axis so it displays the full date
      xlab("Timestamp (D/M/Y H:M:S)") +
      ggtitle(paste0(i,' - ', unique(markerSigs$age))) +
       geom_point(aes(ts, sig, color = antBearing, label = ts, text = runID)) + 
      geom_point(data = markerFlights, aes(ts, sig), size = 4, shape = 1, stroke = 2) + 
      guides(color = FALSE) +
      facet_grid(site+antBearing~.)
      saveWidget(ggplotly(p, tooltip = c("label", "text")), paste0('proj#',projectID,'-tagPlot#', i ,'.html'))
     
    }
}


i <- '1861-32426'
tag_df <- fixed_tagList %>% filter(markerNumber == i) %>% mutate(site = ifelse(is.na(site), ifelse(is.na(recv), 'UNK', recv), site)) %>%
  mutate(ts = as.POSIXct(ts, origin = '1970-01-01'), antBearing = as.factor(antBearing))

## Check if it there is any data
if (nrow(tag_df) > 0 & !is.na(i)) {
  
  sunRiseSetTimes <- tibble(ts = days <- seq(min(tag_df$ts)-days(1), max(tag_df$ts), 24*3600),
                            lat = rep(43.4674, length(ts)),
                            lon = rep(-65.7489, length(ts)))
  
  sunRiseSetTimes <- sunRiseSetTimes %>% sunRiseSet(lat = "lat", lon = "lon", ts = "ts") %>%
    mutate(date = as.Date(ts)) %>%
    select(-lat, -lon, -ts)
  
  attr(sunRiseSetTimes$sunset, 'tzone') <- "America/Halifax"
  attr(sunRiseSetTimes$sunrise, 'tzone') <- "America/Halifax"
  
  bandingDates <- rawBandingData$DATE.BANDED[rawBandingData$markerNumber == i]
  
  if (nrow(tag_df %>% filter(year(ts) == year(min(bandingDates)))) > 0) {
    bandingDates <- bandingDates[year(bandingDates)==min(year(bandingDates))]
    tag_df <- tag_df %>% filter(year(ts) == year(min(bandingDates)))
  } else {
    bandingDates <- bandingDates[year(bandingDates)==max(year(bandingDates))]
    tag_df <- tag_df %>% filter(year(ts) == year(max(bandingDates)))
  }
  
  latestDetection <- ifelse(max(bandingDates) > max(tag_df$ts), max(bandingDates), max(tag_df$ts))
   p <- plot_ly(data = tag_df, x = ~ts, y = ~sig, color = ~antBearing) %>%
    add_trace()
  p
}













## Make a function that will chose an appropriate filename and scale based on number of plots
savePlot <- function (markerNumber) {
  print(markerNumber)
  if (!is.null(p[[markerNumber]]) & !is.na(markerNumber)) {
    numPlots <- length(p[[markerNumber]])
    if (numPlots > 1) toPlot <- do.call(arrangeGrob, p[[markerNumber]]) else toPlot <- p[[markerNumber]][[1]]
    ggsave(paste0('proj#',projectID,'-tagPlot#',markerNumber,'.png'), plot = toPlot, scale = 0.5 + floor(sqrt(numPlots)))
  }
}

############ SAVE PLOTS ##############
## Run this to loop through each plot and save the images!
mapply(savePlot, unique(fixed_tagList$markerNumber))
######################################
