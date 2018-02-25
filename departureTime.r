## Load libraries

library(tidyverse)
library(devtools)
library(RSQLite)
library(grid)
library(gridExtra)

library(rworldmap)
library(ggmap)

########
### Import datasets into list
########
if (!exists('projectID')) projectID <- 109

source("getMetaData.r")
source("loadMotusData.r")

bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')


### Load TAG DETECTIONS ###
tag.sql <- dplyr::src_sqlite(paste0('project-',projectID,'.motus'))
tag.tbl <- tbl(tag.sql, "alltags")
rawData <- hits.tbl %>%  
  distinct %>% collect %>% as.data.frame 

### Load BANDING METADATA ###
rawBandingData <- loadMeta()
colnames(rawBandingData)[2] <- "markerNumber"

### Load RECEIVER METADATA ###
receiverData <- read_csv("data/receiver-deployments.csv") %>%
  mutate(
    tsStart = as.POSIXct(tsStart,origin='1970-01-01'),
    tsEnd = as.POSIXct(tsEnd,origin='1970-01-01'),
    deploymentName = gsub("<[^\\]]*>", "", deploymentName, perl = TRUE)
  )

receiverData[receiverData$receiverID == 'SG-F783RPI36B63',]
# Load Shitlist
shitList <- read_csv('109-shitSites.csv', col_names = c('site'))$site


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
    weightByWing = WEIGHT/WING
  )
bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')

# Fix raw data dates/times and ID
rawDataComb <- rawData %>%
  #  filter(lat != 0, lon != 0) %>%
  #  rename(gpsLat = lat, gpsLon = lon, site = recvDeployName)
  mutate(
    ts = as.POSIXct(ts,origin='1970-01-01'),
    #    id = as.character(motusTagID),
    date = format(ts,"%Y-%m-%d"),
    site = ifelse(site %in% bpsites, 'Bon Portage', ifelse(site %in% sealsites, 'Seal Island', site))
  ) %>%
  left_join(bData, by = 'markerNumber')

# Filter receivers by those available while tags were acive
receiverLocations <- receiverData %>%
  filter(is.na(dtEnd) | tsEnd > min(rawData$ts),
         tsStart < max(rawData$ts)+days(365)) %>%
  # Select relevant columns
  select(recv = recvDeployID, lon = longitude, lat = latitude, site = deploymentName) %>%
  mutate(recv = as.character(recv))

# Create new summary tibble of receivers
allSitesWithDetections <- rawData %>%
  group_by(recv) %>%
  summarise(lat = lat[1], site = site[1]) %>%
  filter(!is.na(lat) & site != 'Old Cut' & site != 'OLDCUT')

fixedData <- rawDataComb %>%
  select(id, motusTagID, markerNumber, recv, ts, date, site,
         species = SPECIES.NAME, age = AGE.ch, sex = SEX.ch, wing = WING, weight = WEIGHT, depLon, depLat, lon, lat, runLen, spEN) %>%
  filter(!is.na(markerNumber), !site %in% c('Old Cut','OLDCUT', shitList), runLen > 2, !lon %in% c(NA & 0)) %>%
  mutate(
    date = format(ts,"%Y-%m-%d"),
    hour = format(ts,"%Y-%m-%d %H"),
    species = ifelse(is.na(species),  ifelse(is.na(spEN), 'Unknown', spEN), species), 
    age = ifelse(is.na(age), 'Unknown', age),
    year = format(ts, '%Y'),
    recvDeployName = site,
    tagDeployID = markerNumber
  )


fixedData <- loadMotusData(109)

sitewise <- fixedData %>% siteTrans(latCoord = ~lat, lonCoord = ~lon)

depLLs <- fixedData %>%
  group_by(markerNumber) %>%
  summarise(depLat = depLat[1], depLon = depLon[1])

tagDeps <- sitewise %>%
  rename(markerNumber = tagDeployID) %>%
  left_join(bData, by = 'markerNumber') %>%
  left_join(depLLs, by = 'markerNumber') %>%
  arrange(ts.x) %>%
  group_by(markerNumber) %>%
  summarise(depDate = last(ts.x), arrDate = last(ts.y), arrSite = last(recvDeployName.y), tot_ts = last(tot_ts),
            dist = last(dist)/1000, age = AGE.ch[1], depLon = depLon[1])

tagDeps %>%
  filter() %>%
  mutate(depYear = year(depDate), 
         depDate = as.Date(format(depDate, '%m-%d'),'%m-%d'),
         depLoc = ifelse(depLon > -66, 'Bon Portage Island', 'Seal Island')) %>%
  ggplot(aes(depDate, dist, color = age))+
  geom_jitter()+
  ggtitle('Bon Portage and Seal Island')+
#  ggtitle('Seal Island')+
  facet_grid(depYear~depLoc)+
  scale_x_date(date_labels='%b-%d', date_breaks = "1 week")+
  scale_y_log10()+
  xlab('Departure Date')+
  ylab('Distance to fist site after departure (Km)')

runList <- fixedData %>%
  group_by(runID) %>%
  summarise(sig = max(sig), ts = as.POSIXct(max(ts), origin = '1970-01-01'))

depTimes <- read_csv("Project109-sigplots-Notes.csv")[1:13] %>%
  filter(is.na(died), is.na(nodata)) %>%
  #select() %>%
  left_join(bData, by = 'markerNumber') %>%
  mutate(endBlip = as.POSIXct(endBlip, format = "%m/%d/%Y", origin = '1970-01-01'),
         BP = as.POSIXct(BP, format = "%m/%d/%Y", origin = '1970-01-01'),
         SEAL = as.POSIXct(SEAL, format = "%m/%d/%Y", origin = '1970-01-01')
  ) %>%
#  gather(BP, SEAL, key = 'site', value = 'ts') %>%
  left_join(runList, by = 'runID') %>%
  filter(!is.na(ts))

fixed_tagList %>%
  filter(hitID %in% depFlights$hitID) 
depFlights %>%  
  mutate(ts = as.numeric(format(as.Date(ts),'%j'))) %>%
  ggplot()+
  geom_jitter(aes(ts, bandsite, color = bandsite), height = 0.2, width = 0)+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b"))+
  ggtitle(paste0('Estimated Departure times from banding site for each individual (n = ',length(unique(depFlights$markerNumber)),')', collapse = ''))+
  facet_grid(age~.)

depFlights %>%
  filter(year(ts) == 2017) %>%
  mutate(ts = as.POSIXct(paste(hour(ts),":",minute(ts), sep = ''), format = '%H:%M', origin = '2018-01-01')) %>%
  ggplot()+
  geom_jitter(aes(ts, bandsite, color = bandsite), height = 0.2, width = 0)+
  scale_x_datetime()+
  ggtitle(paste0('Estimated Departure times from banding site for each individual (n = ',length(unique(depFlights[year(depFlights$ts) == 2017,]$markerNumber)),')', collapse = ''))+
  facet_grid(age~.)


depTimes %>%
  filter(year(ts) == 2017) %>%
  mutate(ts = as.POSIXct(paste(hour(ts),":",minute(ts), sep = ''), format = '%H:%M', origin = '2018-01-01')) %>%
  ggplot()+
  geom_jitter(aes(ts, bandsite, color = bandsite), height = 0.2, width = 0)+
  scale_x_datetime()+
  ggtitle(paste0('Estimated Departure times from banding site for each individual (n = ',length(unique(depFlights[year(depFlights$ts) == 2017,]$markerNumber)),')', collapse = ''))+
  facet_grid(age~.)


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

depTimes <- flights %>% 
  left_join(flight.meta, by = 'markerNumber') %>% 
  filter(code != 'X' | is.na(code)) %>% 
  mutate(isDep = runID %in% depFlights$runID) %>%
  left_join(plyr::count(flights, "markerNumber"), by = "markerNumber")

depTimes 

sunRiseSetTimes <- tibble(ts = days <- seq(min(depTimes$ts)-days(1), max(depTimes$ts), 24*3600),
                          lat = rep(43.4674, round(max(depTimes$ts)-min(depTimes$ts))+1),
                          lon = rep(-65.7489, round(max(depTimes$ts)-min(depTimes$ts))+1))

sunRiseSetTimes <- sunRiseSetTimes %>% sunRiseSet(lat = "lat", lon = "lon", ts = "ts") %>%
  mutate(date = as.Date(ts)) %>%
  select(-lat, -lon, -ts)

attr(sunRiseSetTimes$sunset, 'tzone') <- "America/Halifax"

depBySunset <- depTimes %>%
  mutate(date = as.Date(ts), 
         jdate = as.numeric(format(date, "%j")),
         #sunset = sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) > -23*3600, ]$sunset,
         sunset = (lapply(ts, function(ts){sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) > -24*3600, ]$sunset})),
         #sunset = sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) > -24*3600, ]$sunset,
         bySunset = (as.numeric(ts)/60) - (as.numeric(sunset)/60))

#sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(depBySunset$ts[2])) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(depBySunset$ts[2])) > -23*3600, ]$sunset
#sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(depBySunset$ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(depBySunset$ts)) > -23*3600, ]$sunset
#with(depBySunset, sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600, ]$sunset)
lapply(depBySunset$ts, function(ts){sunRiseSetTimes[(as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) < 3600 & (as.numeric(sunRiseSetTimes$sunset)-as.numeric(ts)) > -23*3600, ]$sunset})

depStats <- depBySunset %>% group_by(markerNumber) %>% summarise(age = age[1], nflights = n()) %>% group_by(age) %>% summarise(a = mean(nflights), sd = sd(nflights), nflights = sum(nflights), n = n())

depBySunset %>% group_by(markerNumber) %>% summarise(age = age[1], nflights = n()) %>% ggplot(aes(age, nflights))+geom_boxplot()

depBySunset %>%
  #  filter(year(ts) == 2016) %>%
  ggplot()+
  geom_jitter(aes(bySunset, sig, color = bandsite), height = 0.2, width = 0)+
  ggtitle(paste0('Departures from banding site relative to sunset times for each individual (n = ',length(unique(depBySunset$markerNumber)),')', collapse = ''))+
  xlab("Minutes from Sunset")+
  ylab("Departure site")+
  facet_grid(age~.)

depBySunset %>%
  filter(isDep) %>%
  #  filter(year(ts) == 2016) %>%
  ggplot()+
  scale_x_continuous()+
  geom_jitter(aes(bySunset, bandsite, color = bandsite), height = 0.2, width = 0)+
  ggtitle(paste0('Departures from banding site relative to sunset times for each individual (n = ',length(unique(depBySunset$markerNumber)),')', collapse = ''))+
  #xlab("Departure Date")+
  ylab("Departure site")+
  facet_grid(age~.)

depBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bandsite, bySunset))+
  geom_boxplot()+
  ylab("Minutes from Sunset")+
  xlab("Departure site")+
  facet_grid(age~.)

depBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bandsite, as.numeric(format(date, "%j"))))+
  geom_boxplot()+
  ylab("Julien date")+
  xlab("Departure site")+
  facet_grid(age~.)

depBySunset %>%
  filter(isDep) %>%
  ggplot(aes(age, as.numeric(format(date, "%j"))))+
  geom_boxplot()+
  ylab("Julien date")+
  xlab("Departure site")+
  facet_grid(bandsite~.)

depBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bySunset))+
  geom_histogram()+
  ylab("# Birds")+
  xlab("Minutes after sunset")+
  facet_grid(age~bandsite)

depBySunset %>%
  filter(isDep) %>%
  ggplot(aes(bySunset))+
  geom_histogram()+
  ylab("# Birds")+
  xlab("Minutes after sunset")+
  facet_grid(bandsite~.)

depBySunset %>%
  #filter(freq == 11, bandsite == 'BP') %>% 
  group_by(markerNumber) %>%
  summarise(freq = n(), age = age[1], bandsite = bandsite[1]) %>%
  ggplot(aes(freq))+
  geom_histogram()+
  facet_grid(age~bandsite)F
  geom_jitter()

####################
####################

library(stats)

par(mfrow=c(2,2))

glm1 <- depBySunset %>% filter(isDep)  %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = gaussian)
plot(glm1)
glm2 <- depBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = Gamma)
plot(glm2)

summary(glm2)

anova(glm2, test = 'F')


glm3 <- depBySunset %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = gaussian)
plot(glm3)
glm4 <- depBySunset %>% glm(formula = log(bySunset) ~ age + jdate + bandsite, family = Gamma)
plot(glm4)

summary(glm3)

anova(glm3, test = 'F')


glm5 <- depBySunset %>% glm(formula = jdate ~ age * bandsite, family = gaussian)
plot(glm5)
summary(glm5)
anova(glm5, test = 'F')


glm6 <- depBySunset %>% filter(isDep) %>% glm(formula = jdate ~ age * bandsite, family = gaussian)
plot(glm6)
summary(glm6)
anova(glm6, test = 'F')

glm7 <- depBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age * bandsite, family = gaussian)
plot(glm7)
summary(glm7)
anova(glm7, test = 'F')




glm1 <- depBySunset %>% glm(formula = bySunset ~ age, family = gaussian)
glm1 <- depBySunset %>% glm(formula = bySunset ~ bandsite, family = gaussian)
glm2 <- depBySunset %>% glm(formula = bySunset ~ age * bandsite, family = gaussian)
glm3 <- depBySunset %>% glm(formula = bySunset ~ age + FAT + bandsite, family = gaussian)
glm4 <- depBySunset %>% glm(formula = bySunset ~ age + date + bandsite, family = gaussian)
glm1 <- depBySunset %>% glm(formula = log(bySunset) ~ age + FAT + date + bandsite, family = gaussian)
glm2 <- depBySunset %>% filter(isDep)  %>% glm(formula = log(bySunset) ~ age + date + bandsite, family = gaussian)
glm2 <- depBySunset %>% filter(isDep) %>% glm(formula = log(bySunset) ~ age + date + bandsite, family = Gamma)
glm3 <- depBySunset %>% glm(formula = log(bySunset) ~ age + bandsite, family = Gamma)
glm4 <- depBySunset %>% filter(!is.na(weightByWing)) %>% glm(formula = log(bySunset) ~ age + weightByWing + bandsite, family = gaussian)
glm3
AIC(glm1, glm2, glm6, glm7)
AIC(glm3, glm4, glm5)



?glm

wilcox.test(depBySunset[depBySunset$AGE.ch == 'Hatch-year' & depBySunset$bandsite != 'BP',]$bySunset, depBySunset[depBySunset$AGE.ch == 'Hatch-year' & depBySunset$bandsite == 'BP',]$bySunset, alternative = 'less')

wilcox.test(depBySunset[depBySunset$AGE.ch != 'Hatch-year' & depBySunset$bandsite == 'BP',]$bySunset, depBySunset[depBySunset$AGE.ch == 'Hatch-year' & depBySunset$bandsite == 'BP',]$bySunset, alternative = 'less')

wilcox.test(depBySunset[depBySunset$bandsite != 'BP',]$bySunset, depBySunset[depBySunset$bandsite == 'BP',]$bySunset, alternative = 'less')

t.test(depBySunset[depBySunset$AGE.ch != 'Hatch-year',]$bySunset, depBySunset[depBySunset$AGE.ch == 'Hatch-year',]$bySunset, alternative = 'less')

wilcox.test(depBySunset[depBySunset$AGE.ch != 'Hatch-year',]$bySunset, depBySunset[depBySunset$AGE.ch == 'Hatch-year',]$bySunset, alternative = 'less')

var(depBySunset[depBySunset$AGE.ch == 'Hatch-year' & depBySunset$bandsite != 'BP',]$bySunset)

sitewise %>% filter(site == 'Baccaro') %>% group_by(markerNumber) %>% summarise()
sitewise %>% group_by(site) %>% summarise()

depTimes %>%
  mutate(ts = as.POSIXct(paste(hour(ts),":",minute(ts), sep = ''),format = '%H:%M', origin = '2018-01-01')) %>%
  select(ts)

format(depTimes$ts[1], "%H:%M")

deps <- depTimes %>%
  select(markerNumber, depSite = bandsite, age = AGE.ch)

sitewise <- fixedData %>% 
  group_by(runID) %>%
  summarise(site = site[1], 
            markerNumber = markerNumber[1], 
            lat = lat[1],
            lon = lon[1], 
            ts = as.POSIXct(max(ts), origin='1970-01-01'), 
            motusTagID = motusTagID[1]) %>%
  left_join(bData, by = 'markerNumber') %>%
  mutate(AGE.ch = ifelse(is.na(AGE.ch), 'UNK', ifelse(AGE.ch == 'Hatch-year', AGE.ch, 'Adult')))

depTimes2 <- sitewise %>% 
  filter(site %in% c('Seal Island', 'Bon Portage')) %>%
  left_join(deps, by = 'markerNumber') %>%
#  filter(depSite == 'SEAL' & site == 'Seal Island') %>%
  filter(depSite == 'BP' & site == 'Bon Portage') %>%
 # filter(site == 'BP') %>%
  group_by(markerNumber) %>%
  summarise(depTime = max(ts), age = AGE.ch[1], depSite = depSite[1]) %>%
  arrange(depTime)

sunRiseSetTimes <- tibble(ts = days <- seq(min(depTimes2$depTime), max(depTimes2$depTime), 24*3600),
                          lat = rep(43.4674, round(max(depTimes2$depTime)-min(depTimes2$depTime))),
                          lon = rep(-65.7489, round(max(depTimes2$depTime)-min(depTimes2$depTime)))) %>%
  mutate(ts = as.POSIXct(format(ts, '%d-%m %H:%M'), format = '%d-%m %H:%M', origin = '2018-01-01'))

sunRiseSetTimes <- sunRiseSetTimes %>% sunRiseSet(lat = "lat", lon = "lon", ts = "ts")

depTimes2 %>% 
#  mutate(depTime = as.numeric(format(as.Date(depTime),'%j'))) %>%
  mutate(depTime = as.POSIXct(strftime(depTime,'%d-%m %H:%M'), format = '%d-%m %H:%M')) %>%
  filter(depTime < as.POSIXct('2018-09-17', format = '%Y-%m-%d')) %>%
  ggplot()+
  geom_jitter(aes(depTime, age), height = 0.2, width = 0)+
  scale_x_datetime()+
#  scale_x_continuous(labels = function(x) format(as.Date(as.character(x)), "%d-%b"))
 # scale_x_continuous(labels = function(x) format(x, "%d-%b"))
#  geom_polygon(data = sunRiseSetTimes, aes(x = sunrise), color = '#FF0000')+
  geom_vline(data = sunRiseSetTimes, aes(xintercept = sunrise), color = '#FF0000')+
  geom_vline(data = sunRiseSetTimes, aes(xintercept = sunset), color = '#0000FF', stroke = 'dashed')

names(tagList)

tagList %>% filter(runID == 18628060)
temp <- tagList %>%
#  filter(markerNumber == '1861-32779', site %in% bpsites) %>%
  filter(markerNumber == '1861-32243', site %in% sealsites) %>%
  group_by(runID, ant, site) %>%
  summarise(sig = max(sig), ts = as.POSIXct(max(ts), origin = '1970-01-01')) #%>%
#  filter(year(ts) >= 2017)
temp %>%
  ggplot(aes(ts, sig))+
  geom_jitter()
temp %>%
  ggplot(aes(ts, sig, color = as.factor(ant)))+
  geom_jitter()+
  facet_grid(site~.)

tagList %>%
  filter(markerNumber == '1861-32779', site %in% bpsites) %>%
  mutate(ts = as.POSIXct(ts, origin = '1970-01-01')) %>%
  filter(day(ts) > 1 & month(ts) < 9) %>%
  ggplot(aes(ts, sig, color = as.factor(ant)))+
  geom_jitter()+
  facet_grid(site~.)

fixed_tagList$ts


