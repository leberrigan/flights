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

source("loadMotusData.r")

########
### Import datasets into list
########

# Set timezone to GMT
Sys.setenv(TZ = 'GMT')

# Define Project ID
projectID <- 109

# Local address of database
database <- "data/"

# Create a new SQLITE database?
newDatabase <- F

# Load tag hits
tagHits <- loadMotusData(projectID, database, newDatabase)

# Create a vectors of BP and Seal sites
bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')

# Only select receivers from BP and Seal
localHits <- tagHits %>%
  filter(site %in% c(bpsites, sealsites))

# Load list of flight runIDs
flights.raw <- read_csv("Project109-sigplots-Flights.csv") 

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
         postDep = (runID > depRunID)) %>%
  left_join(plyr::count(flights, "markerNumber"), by = "markerNumber")

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

# Get some stats on the number of flights
flightStats <- flightBySunset %>%
  group_by(markerNumber) %>% 
  summarise(age = age[1], nflights = n()) %>% 
  group_by(age) %>% 
  summarise(a = mean(nflights), sd = sd(nflights), nflights = sum(nflights), n = n())

flightStats

####
#  All flights including departures, not including post-departure flights
####

flightBySunset %>% 
  group_by(markerNumber) %>% 
  summarise(age = age[1], nflights = n()) %>% 
  ggplot(aes(age, nflights))+
  geom_boxplot()

flightBySunset %>%
  ggplot()+
  geom_jitter(aes(bySunset, sig, color = bandsite), height = 0.2, width = 0)+
  ggtitle(paste0('Departures from banding site relative to sunset times for each individual (n = ',length(unique(flightBySunset$markerNumber)),')', collapse = ''))+
  xlab("Minutes after Sunset")+
  ylab("Departure site")+
  facet_grid(age~.)

flightBySunset %>%
  group_by(markerNumber) %>%
  summarise(freq = n(), age = age[1], bandsite = bandsite[1]) %>%
  ggplot(aes(freq))+
  geom_histogram()+
  facet_grid(age~bandsite)+
  xlab("Number of flights")+
  ylab("Number of birds")


####
#  All departures
####

flightBySunset %>%
  filter(isDep) %>%
  ggplot()+
  scale_x_continuous()+
  geom_jitter(aes(bySunset, bandsite, color = bandsite), height = 0.2, width = 0)+
  ggtitle(paste0('Departures from banding site relative to sunset times for each individual (n = ',length(unique(flightBySunset$markerNumber)),')', collapse = ''))+
  xlab("Minutes after sunset")+
  ylab("Departure site")+
  facet_grid(age~.)

flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bandsite, bySunset))+
  geom_boxplot()+
  ylab("Minutes after Sunset")+
  xlab("Departure site")+
  facet_grid(age~.)


flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bandsite, as.numeric(format(date, "%j"))))+
  geom_boxplot()+
  ylab("Julien date")+
  xlab("Departure site")+
  facet_grid(age~.)

flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(age, as.numeric(format(date, "%j"))))+
  geom_boxplot()+
  ylab("Julien date")+
  xlab("Departure site")+
  facet_grid(bandsite~.)

flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bySunset))+
  geom_histogram()+
  ylab("# Birds")+
  xlab("Minutes after sunset")+
  facet_grid(age~bandsite)

flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bySunset))+
  geom_histogram()+
  ylab("# Birds")+
  xlab("Minutes after sunset")+
  facet_grid(bandsite~.)

  
####################
# Models
####################

# Change the viewport so that it can display 4 plots at once
par(mfrow=c(2,2))

glm1 <- flightBySunset %>% filter(isDep)  %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = gaussian)
plot(glm1)
glm2 <- flightBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = Gamma)
plot(glm2)

summary(glm2)

anova(glm2, test = 'F')


glm3 <- flightBySunset %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = gaussian)
plot(glm3)
glm4 <- flightBySunset %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = Gamma)
plot(glm4)

summary(glm3)

anova(glm3, test = 'F')


glm5 <- flightBySunset %>% glm(formula = jdate ~ age * bandsite, family = gaussian)
plot(glm5)
summary(glm5)
anova(glm5, test = 'F')


glm6 <- flightBySunset %>% filter(isDep) %>% glm(formula = jdate ~ age * bandsite, family = gaussian)
plot(glm6)
summary(glm6)
anova(glm6, test = 'F')

glm7 <- flightBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age * bandsite, family = gaussian)
plot(glm7)
summary(glm7)
anova(glm7, test = 'F')




glm1 <- flightBySunset %>% glm(formula = bySunset ~ age, family = gaussian)
glm1 <- flightBySunset %>% glm(formula = bySunset ~ bandsite, family = gaussian)
glm2 <- flightBySunset %>% glm(formula = bySunset ~ age * bandsite, family = gaussian)
glm3 <- flightBySunset %>% glm(formula = bySunset ~ age + FAT + bandsite, family = gaussian)
glm4 <- flightBySunset %>% glm(formula = bySunset ~ age + date + bandsite, family = gaussian)
glm1 <- flightBySunset %>% glm(formula = log(bySunset) ~ age + FAT + date + bandsite, family = gaussian)
glm2 <- flightBySunset %>% filter(isDep)  %>% glm(formula = log(bySunset) ~ age + date + bandsite, family = gaussian)
glm2 <- flightBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age + date + bandsite, family = Gamma)
glm3 <- flightBySunset %>% glm(formula = log(bySunset) ~ age + bandsite, family = Gamma)
glm4 <- flightBySunset %>% filter(!is.na(weightByWing)) %>% glm(formula = log(bySunset) ~ age + weightByWing + bandsite, family = gaussian)
glm3
AIC(glm1, glm2, glm6, glm7)
AIC(glm3, glm4, glm5)

