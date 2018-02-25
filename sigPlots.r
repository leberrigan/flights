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

### Load TAG DETECTIONS ###
t <- dplyr::src_sqlite(paste0('project-',projectID,'.motus'))
df <- tbl(t, "alltags")

# Load Shitlist
shitList <- read_csv('109-shitSites.csv', col_names = c('site'))$site


# Collect all tag observations and flatten them into a tibble
tagList <- df %>% collect %>% as_tibble %>% filter(runLen > 2)


### Load BANDING METADATA ###
source("getMetaData.r")
rawBandingData <- loadMeta()
colnames(rawBandingData)[2] <- "markerNumber"

bData <- rawBandingData %>% 
  group_by(markerNumber) %>%
  arrange(DATE.BANDED) %>%
  summarise(
    SPECIES.NAME = SPECIES.NAME[1],
    SPECIES.CODE = SPECIES.CODE[1],
    AGE.ch = checkAge(AGE.ch),
    SEX.ch = checkSex(SEX.ch),
    WING = mean(WING),
    FAT = (max(FAT)),
    WEIGHT = mean(WEIGHT),
    TIMES.TRAPPED = n(),
    AGE = AGE[1],
    weightByWing = WEIGHT/WING,
    DATE.BANDED = min(DATE.BANDED)
  )

# If the 'id' column is named 'mfgID', rename it
if ('mfgID' %in% colnames(tagList)) {
  tagList <- tagList %>%
  mutate(id = ifelse(is.na(mfgID), motusTagID, mfgID))
}

## Load the receiver deployments and then rbind tibbles
# You will have to edit this part depending on what deployment file you are using
# Summarise by receiver ID and select relevant columns based on latest deployment
receivers <- read_csv('data/receiver-deployments.csv') %>% #select(-motusDeviceID, -test, -fixtureType)%>%
  group_by(receiverID) %>%
  summarise(latitude = latitude[recvDeployID == max(recvDeployID)][1], longitude = longitude[recvDeployID == max(recvDeployID)][1], deploymentName = deploymentName[recvDeployID == max(recvDeployID)][1]) %>%
  rename(recv = receiverID)

bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')

#unique(tagList$recv) == "SG-4F30RPI36C9C"
#test <- filter(fixed_tagList, recv == "SG-4F30RPI36C9C", markerNumber == '1861-32705')




source("loadMotusData.r")
tagList <- loadMotusData(109)
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


test<- tibble(mn = unique(fixed_tagList$markerNumber), died = NA, endBlip = NA)
#write_csv(test, 'markerNotes.csv')

######################################
## Data Validation
######################################

# Make a list of sites that are known to be false
suspectSites <- shitSites
  c('Canopy Tower', 'La Victoria 1', 'La Victoria 2', 
                  'LLICALDAD', 'Los Vientos Forest', 'Quempillen (Chile)',
                  'Koffler', 'Old Cut', 'GB3', 
             #     'BPhill', 'BPLH', 'BPnorthomni', 'BPwestomni', 
                  'Grand-Ile')

## Group tag list by runID and get maximum burst slop for each run
# runID is the ID of the particular sequence tag pulses
# Burst slop is the difference between the measured burst rate and the actual tag burst rate
tagRuns <- fixed_tagList %>%
  arrange(ts) %>% 
  group_by(id, site, runID) %>% 
  summarise(ant = ant[1], 
            tsStart = min(ts), # The start of the run
            runLen = runLen[1], 
            burstSlopMax = max(abs(burstSlop)), 
            freqsdMax = max(freqsd),
            lat = ifelse(median(lat)!=0, median(lat), max(lat)),
            lon = ifelse(median(lon)!=0, median(lon), min(lon))) %>%
  mutate(suspect = site %in% suspectSites) # Suspected sites will be TRUE

# Plot the burstSlopMax over tsStart and grid it by site (lump tags together)
ggplot(tagRuns, aes(color = as.factor(ant))) +
  geom_point(aes(as.POSIXct(tsStart, origin = '1970-01-01'), burstSlopMax, color = burstSlopMax))+
  scale_x_datetime(date_labels = "%d/%m/%y %H:%M:%S") + # Format the x-axis so it displays the full date
  xlab("Timestamp (D/M/Y H:M:S)") +
  ggtitle(paste0(with(tagRuns,site), ' Tag 44.1 Problem Data')) +
  guides(color = FALSE) +
  facet_grid(.~site)

# Collapse into sites
# Create new variable 'Valitity' based on freqsd and suspect 
tagRuns_bysite <- tagRuns %>%
  group_by(site) %>%
  summarise(burstSlopMax = median(burstSlopMax),
            suspect = suspect[1], 
            freqsdMax = median(freqsdMax), 
            numAnt = length(unique(ant)),
            lat = lat[1],
            lon = lon[1],
            tsStart = tsStart[1]) %>%
  mutate(validity = ifelse(suspect, ifelse(freqsdMax < 0.1, 'Suspect + freqsd < 0.1', 'Suspect + freqsd >= 0.1'), ifelse(freqsdMax < 0.1, 'Not suspect + freqsd < 0.1', 'Not suspect + freqsd >= 0.1')))


# Plot the burstSlopMax over site and color it based on validity (lump tags together)
ggplot(tagRuns_bysite, aes(site, burstSlopMax, color = as.factor(validity), shape = as.factor(numAnt))) +
  guides(color = guide_legend(title = NULL))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(size = 3)

# Make a boxplot to show sites where burstSlopMax > 0.01
ggplot(filter(tagRuns, burstSlopMax > 0.01), aes(site, burstSlopMax, color = suspect)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Make a new high resolution map
worldMap <- getMap(resolution = "high")
# Connect up all the points so the polygons are closed
worldMap.df <- fortify(worldMap)
worldMap <- NULL

# Create a new suspect list of sites where burstSlopMax > 0.01
newSuspectList <- tagRuns_bysite[tagRuns_bysite$burstSlopMax>0.01, ]

# Get lat/lon bounding box around these sites
latLonBounds <- with(newSuspectList[!is.na(newSuspectList$lat),], list(
                                        c((mean(lon)-(diff(range(lat))/4)),(mean(lon)+(diff(range(lat))/4))),
                                        range(lat)
                                      )
                    )

# Map new suspect data
ggplot(newSuspectList[!is.na(newSuspectList$lat),], aes(lon, lat, fill = as.factor(validity)))+
  geom_polygon(data = worldMap.df, aes(long, lat,group=group), fill="#444444", colour="#999999")+
  geom_point(# Plot all receivers as empty black circles
    shape = 21,
    stroke = 1, 
    size = 3, 
    alpha = 0.7)+
  coord_fixed(xlim = latLonBounds[[1]], ylim = latLonBounds[[2]])+
  theme(panel.background = element_rect(fill = '#222222'),
        panel.grid.major = element_line(colour = '#666666', linetype = 'dashed'),
        panel.grid.minor = element_line(colour = '#666666', linetype = 'dotted'),
        legend.key = element_rect(fill = '#888888'))


# Create a new good data list of sites where burstSlopMax <= 0.01
goodList <- tagRuns_bysite[tagRuns_bysite$burstSlopMax<=0.01, ]

# Get lat/lon bounding box around these sites
latLonBounds <- with(goodList[!is.na(goodList$lat),], list(
                                                            range(lon),
                                                            range(lat)
                                                          )
                    )

# Map new good data
ggplot(goodList[!is.na(goodList$lat),], aes(lon, lat, fill = as.factor(validity)))+
  geom_polygon(data = worldMap.df, aes(long, lat,group=group), fill="#444444", colour="#999999")+
  geom_point(# Plot all receivers as empty black circles
    shape = 21,
    stroke = 1, 
    size = 3, 
    alpha = 0.7)+
  coord_fixed(xlim = latLonBounds[[1]], ylim = latLonBounds[[2]])+
  theme(panel.background = element_rect(fill = '#222222'),
        panel.grid.major = element_line(colour = '#666666', linetype = 'dashed'),
        panel.grid.minor = element_line(colour = '#666666', linetype = 'dotted'),
        legend.key = element_rect(fill = '#888888'))
