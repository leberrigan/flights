library(tidyverse) 
library(motus)


loadMotusData <- function (projectID) {

hits.sql <- tagme(projectID, update= TRUE, forceMeta=TRUE)

hits.tbl <- tbl(hits.sql, "alltags")

tagMeta <- metadata(hits.sql, projectIDs = projectID)

hits.df <- select(hits.tbl, 
                  motusTagID, id, hitID, runID, batchID, 
                  ts, sig, runLen, freqsd, sigsd, slop, burstSlop,
                  antType, antBearing, lat, lon, recv,
                  depLat, depLon, site, markerNumber, spEN) %>% 
  distinct() %>% collect() %>% as.data.frame() %>%
  mutate(ts = as.POSIXct(ts, origin="1970-01-01"), 
         year = year(ts))

### Load RECEIVER METADATA ###
receiverData <- read_csv("data/receiver-deployments.csv") %>%
  mutate(
    tsStart = as.POSIXct(tsStart,origin='1970-01-01'),
    tsEnd = as.POSIXct(tsEnd,origin='1970-01-01'),
    deploymentName = gsub("<[^\\]]*>", "", deploymentName, perl = TRUE)
  ) %>%
  select(recv = receiverID, latitude, longitude, recvDepStart = tsStart, recvDepEnd = tsEnd) 

## Load Shitlist
shitList <- read_csv('109-shitSites.csv', col_names = c('site'))$site

### Load Tag deployment METADATA ###
tagDeploymentData <- read_csv("data/tag-deployments.csv") %>%
  mutate(
    id = as.factor(mfgID),
    tsStart = as.POSIXct(tsStart,origin='1970-01-01'),
    tsEnd = as.POSIXct(tsEnd,origin='1970-01-01'),
    age = ifelse(age == 2, 'Hatch-year', 'Adult')
  ) %>%
  select(markerNumber, tsStart, tsEnd, age, sex, tagDeployID, weight, wing) 

# Fix raw data dates/times and ID
rawDataComb <- hits.df %>%
  #  filter(lat != 0, lon != 0) %>%
  #  rename(gpsLat = lat, gpsLon = lon, site = recvDeployName)
  mutate(
    ts = as.POSIXct(ts,origin='1970-01-01'),
    #    id = as.character(motusTagID),
    date = format(ts,"%Y-%m-%d"),
  #  site = ifelse(site %in% bpsites, 'Bon Portage', ifelse(site %in% sealsites, 'Seal Island', site)),
    lat = ifelse(lat == 0 | is.na(lat), receiverData[receiverData$recv == recv & receiverData$recvDepStart <= ts & (is.na(receiverData$recvDepEnd) | receiverData$recvDepEnd >= ts), ]$latitude,lat),
    lon = ifelse(lat == 0 | is.na(lon), receiverData[receiverData$recv == recv & receiverData$recvDepStart <= ts & (is.na(receiverData$recvDepEnd) | receiverData$recvDepEnd >= ts), ]$longitude,lon)
#    lon = ifelse(lon == 0 | is.na(lon), ,lon)
  ) %>%
  left_join(tagDeploymentData, by = 'markerNumber')



# Fix the data Again
rawDataComb %>%
  select(id,
         hitID,
         runID,
         sig,
         motusTagID, 
         markerNumber, batchID, 
         freqsd, sigsd, slop, burstSlop,
         antType, antBearing, 
         recv,tagDeployID,
         ts, 
         date, 
         site,
         age,
         sex,
         wing,
         weight, 
         depLon, 
         depLat, 
         lon, 
         lat, 
         runLen, 
         spEN) %>%
#  filter(!is.na(markerNumber), !site %in% c('Old Cut','OLDCUT', shitList), runLen > 2, lon != 0 & !is.na(lon)) %>%
  mutate(
    date = format(ts,"%Y-%m-%d"),
    hour = format(ts,"%Y-%m-%d %H"),
    age = ifelse(is.na(age), 'Unknown', age),
    year = format(ts, '%Y'),
    recvDeployName = site,
    tagDeployID = markerNumber
  )

}
