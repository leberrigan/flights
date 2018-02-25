##############################
## Designed for MotusData.r ##
##############################

loadMeta <- function() {
  
  dataDir <- "data/banding/"
  
  fileNames <- c("2016 ABO-FALL BANDED Bandit Aux.CSV",
                "2016 ABO-FALL RETRAPS Bandit Aux.CSV",
                "2016 SUMMER BANDED Bandit Aux.CSV",
                 "2016 SUMMER RETRAPS Bandit Aux.CSV",
                 "ABO FALL 2017 BANDED bandit_FIX.CSV",
                 "ABO FALL 2017 RETRAPS bandit_FIX.CSV",
                 "BANDING-SUMMER-2017-BP bandit_FIX.CSV")
  
  # Load the files and lump them into a single tibble
  splitRawData <- lapply(paste(dataDir, fileNames, sep = ''), function(x) {
    read_csv(x, col_types = cols(DISPOSITION = col_character())) %>%
      separate(TIME.CAUGHT, into = c('HOUR','MINUTE','SECOND')) %>%
      select(DISPOSITION, 
             BAND.NUMBER, 
             SPECIES.NAME,
             SPECIES.CODE, 
             STATUS,
             AGE, 
             SEX, 
             WING, 
             FAT, 
             WEIGHT, 
             DAY, 
             MONTH, 
             YEAR,
             HOUR,
             MINUTE,
             LOCATION)

  })
  
  # Read in IPSP data
  IPSPData <- read_tsv(paste0(dataDir,"IPSP SPRING 2016 BANDED.txt"),
                       col_types = cols(DISPOSITION = col_character(),
                                        DAY = col_integer(),
                                        MONTH = col_integer(),
                                        YEAR = col_integer(),
                                        'HOUR TRAPPED' = col_integer(),
                                       'MINUTE TRAPPED' = col_integer()))
  
  names(IPSPData) <- gsub(" ",".",names(IPSPData),fixed=TRUE)
  
  IPSPData <- IPSPData %>% 
    mutate(
      TIME.CAUGHT = hm(paste(HOUR.TRAPPED, MINUTE.TRAPPED, sep = ':'))
    ) %>%
    select(
      DISPOSITION, 
      BAND.NUMBER, 
      SPECIES.NAME,
      SPECIES.CODE, 
      STATUS,
      AGE, 
      SEX, 
      WING, 
      FAT, 
      WEIGHT, 
      DAY, 
      MONTH, 
      YEAR, 
      TIME.CAUGHT,
      LOCATION
  ) %>% mutate(All_Marker_info_at_release = '')
  
  # Bind all data together  
  rawData <- bind_rows(splitRawData)  
  
  # Replace spaces in variable names with '.'
  names(rawData) <- gsub(" ",".",names(rawData), fixed=TRUE)
  # Create a date object
  rawData <- mutate(rawData, DATE.BANDED = as.POSIXct(paste(DAY, MONTH, YEAR,HOUR, MINUTE, sep = '-'), format="%d-%m-%Y-%H-%M"))
  
  # Codes to Join AGE and SEX
  sexes <- tibble(
    SEX = c(0:5),
    SEX.ch = c('Unknown','','','','Male','Female')
  )
  ages <- tibble(
    AGE = c(0:6),
    AGE.ch = c('Unknown','Adult','Hatch-year','Juvenile','Local','Second-year','After Second-year')
  )
  
  # Prepare SWTH data
  rawData %>%
    left_join(sexes,by = "SEX") %>%
    left_join(ages,by = "AGE")
}

checkSex <- function(sex) {
  matches <- unique(sex[grep("Male|Female", sex)])
  ifelse(length(matches) > 0, ifelse(length(matches) == 1, matches, 'Male?Female'), 'Unknown')
}

checkAge <- function(age) {
  matches <- unique(age[grep("Adult|Hatch-year|Juvenile|Local|Second-year|After Second-year|After Third-year", age)])
  ifelse(length(matches) > 0, matches, 'Unknown')
}

getBandingMeta <- function(bn, bData) { # bn = band number
#  if (nrow(bData[bData$BAND.NUMBER == bn,]) >= 0) {
  bData %>% 
    filter(BAND.NUMBER == bn) %>% 
    group_by(BAND.NUMBER) %>%
    arrange(DATE.BANDED) %>%
    summarise(
      SPECIES.NAME = SPECIES.NAME[1],
      SPECIES.CODE = SPECIES.CODE[1],
      AGE.ch = checkAge(AGE.ch),
      SEX.ch = checkSex(SEX.ch),
      WING = mean(WING),
      FAT = paste(FAT, collapse = ','),
      WEIGHT = mean(WEIGHT),
      TIMES.TRAPPED = n()
    )
}

