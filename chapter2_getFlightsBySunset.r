
#####
##  DEPENDENCIES
#####
# Scripts
#   > loadMotusData.r
# Datasets 
subfolder <- "data/" # default
#   > receiver-deployments.csv
#   > tag-deployments.csv
#   > Project109-sigplots-Flights.csv
#####


## Load libraries

library(tidyverse)
library(devtools)
library(RSQLite)

library(motus)

library(grid)
library(gridExtra)

library(rworldmap)
library(ggmap)

library(stats)
library(lmtest)

source("loadMotusData.r")

########
### Import datasets into list
########

# Set timezone to GMT
Sys.setenv(TZ = 'GMT')

# Define Project ID
projectID <- 109

# Create a new SQLITE database?
newDatabase <- F

# Load tag hits
tagHits <- loadMotusData(projectID, subfolder, newDatabase)

# Create a vectors of BP and Seal sites
bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')

# Only select receivers from BP and Seal
localHits <- tagHits %>%
  filter(site %in% c(bpsites, sealsites))

# Load list of flight runIDs
flights.raw <- read_csv(paste0(subfolder, "Project109-sigplots-Flights.csv", collapse = '')) 

# Get runID for all departure flights
flights.raw <- flights.raw %>%
  mutate(depRunID = ifelse(is.na(depRunID), max(flights.raw[flights.raw$markerNumber == markerNumber, 5:16]), depRunID))

# Select all hits which correspond to each flight runID
flights.df <- localHits[localHits$runID %in% as.numeric(unlist(flights.raw[5:16])),]

# Group flights by runID and select hitID and ts that corresponds to max signals strength during that run
flights <- flights.df %>%
  group_by(runID) %>%
  summarise(ts = ts[which.max(sig)], hitID = hitID[which.max(sig)],sig = max(sig), markerNumber = markerNumber[1], age = age[1])

# Create a new tibble of just departure flights and join it with flight meta data
# This is because some flights were detected after the departure flight
depFlights <- flights %>%
  left_join(flight.meta, by = 'markerNumber') %>%
  group_by(markerNumber) %>%
  summarise(depRunID = ifelse(is.na(depRunID), runID[which.max(ts)], max(depRunID))[1], 
            ts = ts[runID==depRunID], hitID = hitID[runID==depRunID], sig = sig[runID==depRunID], runID = runID[runID==depRunID], 
            bandsite = bandsite[1], age = age[1], code = code[1]) %>%
  select(-depRunID)

# Create new tibble of flight metadata
flight.meta <- depFlights %>% select(bandsite, markerNumber, code, depRunID = runID)

# Create a new tibble of flight times with isDep (T/F), postDep (T/F), and freq (total # flights)
flightTimes <- flights %>% 
  left_join(flight.meta, by = 'markerNumber') %>% 
  filter(code != 'X' | is.na(code)) %>% 
  mutate(isDep = runID %in% depFlights$runID,
         postDep = (runID > depRunID))

# Remove flights that occur after departure
flightTimes <- flightTimes %>% filter(!postDep)

# Create a sequence of dates from the earliest to latest detection date
# Add lat/lon from BP/Seal area
sunRiseSetTimes <- tibble(ts = days <- seq(min(flightTimes$ts)-days(1), max(flightTimes$ts), 24*3600),
                          lat = rep(43.4674, round(max(flightTimes$ts)-min(flightTimes$ts))+1),
                          lon = rep(-65.7489, round(max(flightTimes$ts)-min(flightTimes$ts))+1))

# Calculate sunrise/set times for each date
sunRiseSetTimes <- sunRiseSetTimes %>% sunRiseSet(lat = "lat", lon = "lon", ts = "ts") %>%
  mutate(date = as.Date(ts)) %>%
  select(-lat, -lon, -ts)

# Calculate the flight times as minutes from sunset
flightBySunset <- flightTimes %>%
  mutate(date = as.Date(ts), 
         jdate = as.numeric(format(date, "%j")),
         sunset = (lapply(ts, function(ts){sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) > -24*3600, ]$sunset})),
         bySunset = (as.numeric(ts)/60) - (as.numeric(sunset)/60))

write_csv(flightBySunset %>% select(-sunset), "flightBySunset.csv")