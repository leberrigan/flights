##
##  DEPENDENCIES
##
#
# > flightsBySunset.csv
#
##

## Load Libraries
library(tidyverse)
library(motus)
library(stats)
library(lmtest)
library(survival)

# Set viewport to 2 x 2 grid
par(mfrow=c(2,2))

# Read in flights dataframe
flightBySunset <- read_csv("flightBySunset.csv")

# Get some stats on the number of flights
flightStats <- flightBySunset %>%
  group_by(markerNumber) %>% 
  summarise(age = age[1], nflights = n()) %>% 
  group_by(age) %>% 
  summarise(a = mean(nflights), sd = sd(nflights), nflights = sum(nflights), n = n())

flightStats

flightBySunset %>%
  group_by(markerNumber) %>% 
  summarise(bandsite = bandsite[1], nflights = n()) %>% 
  group_by(bandsite) %>% 
  summarise(a = mean(nflights), sd = sd(nflights), nflights = sum(nflights), n = n())


####
#  MODELS
####

######
### NUMBER OF FLIGHT
######

nFlights <- flightBySunset %>% 
#  filter(!isDep) %>%
  group_by(markerNumber) %>% 
  summarise(age = age[1], bandsite = bandsite[1], n = n()-1)

nFlights %>%
  ggplot(aes(bandsite, n))+
  geom_boxplot()+
  facet_grid(age~.)

nFlights %>%
  ggplot(aes(n)) +
  geom_histogram()+
  facet_grid(age~bandsite)

nFlights %>%
  ggplot(aes(n)) +
  geom_histogram()

sum(nFlights$n)

nFlights %>%
  group_by(age) %>%
  summarise(a = mean(n), sd = sd(n), nInd = n(), n = sum(n))

nFlights %>%
  group_by(bandsite) %>%
  summarise(a = mean(n), sd = sd(n), nInd = n(), n = sum(n))

flightBySunset %>% 
  group_by(markerNumber) %>%
  summarise(age = age[1], bandsite = bandsite[1]) %>%
  group_by(age) %>%
  tally()

flightBySunset %>% 
  group_by(markerNumber) %>%
  summarise(age = age[1], bandsite = bandsite[1], isDep = length(which(isDep==F))) %>% 
  group_by(bandsite) %>% 
  summarise(p = length(which(isDep>0))/n())

nFlights %>%
  group_by(bandsite) %>%
  summarise(n = n())
preDep.flightBySunset <- flightBySunset %>% filter(!isDep)

# Test: Is there a difference in the likelihood of individuals making flights among ages and bandsites?

glm1 <- glm(data = nFlights, formula = (n > 0) ~ age * bandsite, family = binomial)
glm2 <- glm(data = nFlights, formula = (n > 0) ~ age + bandsite, family = binomial)
anova(glm1, glm2, test="Chi")

glm1.1 <- nFlights %>% glm(formula = (n > 0) ~ bandsite, family = binomial)
glm1.2 <- nFlights %>% glm(formula = (n > 0) ~ age, family = binomial)

plot(glm1)

summary(glm1)
anova(glm1, test="Chi")

anova(glm1, update(glm1,.~ -age:bandsite), test = 'Chi')

anova(glm1.2, glm1, test = 'LRT')
anova(glm1.2, glm1, test = 'LRT')

lrtest(glm1.1, glm1)
anova(glm1.1, glm1, test = 'LRT')


# Test: Do individuals at remote sites make more nocturnal flights and
#       do adults make fewer flights than hatch-years and show no difference between sites? 

glm2 <- nFlights %>% glm(formula = n ~ age * bandsite, family = poisson)

plot(glm2)

summary(glm2)

anova(glm2, test = 'F')

lrtest(glm2, glm1)

######
### PRE-DEPARTURE FLIGHT TIMES 
######

flightBySunset %>%
  filter(!isDep) %>% 
  ggplot(aes(bandsite, (bySunset)))+
  geom_boxplot()+
  facet_grid(age~.)+
  ylab('Minutes since sunset')+
  xlab('Banding site')

flightBySunset %>%
  filter(!isDep) %>% 
  ggplot(aes(log(bySunset))) +
  geom_histogram()+
  facet_grid(age~bandsite)

glm3 <- flightBySunset %>% filter(!isDep) %>% glm(formula = log(bySunset) ~ age * bandsite, family = gaussian)

plot(glm3)

summary(glm3)

anova(glm3, test = 'F')

######
### DEPARTURE TIMES (LOG)
######

dFlights <- flightBySunset %>% filter(isDep) %>% select(markerNumber, bySunset, age, bandsite, ts)
  
dFlights %>%
  ggplot(aes(bandsite, (bySunset)))+
  geom_boxplot()+
  facet_grid(age~.)+
  ylab('Minutes since sunset')+
  xlab('Banding site')

dFlights %>%
  ggplot(aes(log(bySunset))) +
  geom_histogram()+
  facet_grid(age~bandsite)

dFlights %>%
  group_by(age) %>%
  summarise(a = mean(bySunset), sd = sd(bySunset), n = n())

dFlights %>%
  group_by(bandsite) %>%
  summarise(a = mean(bySunset), sd = sd(bySunset), n = n())

# Test:  Do individuals at remote sites depart earlier in the evening and
#        do adults depart later than hatch-years and show no difference between sites? 

glm6 <- dFlights %>% glm(formula = log(bySunset) ~ age * bandsite, family = gaussian)

plot(glm6)

summary(glm6)

anova(glm6, test = "F")


######
### NUMBER FLIGHT PRIOR TO DEPARTURE 
######

dDep <- flightBySunset %>%
  rowwise() %>%
  #  mutate(depDate = dFlights[dFlights$markerNumber == markerNumber,]$ts) 
  mutate(daysToDep  = difftime(dFlights[dFlights$markerNumber == markerNumber,]$ts, ts, units = 'days')) 

dDep %>%
  filter(daysToDep > 0) %>%
  ggplot(aes(as.integer(daysToDep), group = age, color= age)) +
  geom_density() +
  facet_grid(.~bandsite)


######
### DEPARTURE DATE
######

flightBySunset %>%
  filter(isDep) %>%
  ggplot(aes(jdate, group = age, color = age)) +
  geom_density() +
  facet_grid(.~bandsite)

glm7 <- flightBySunset %>% filter(isDep) %>% glm(formula = (jdate) ~ age * bandsite, family = gaussian)

plot(glm7)

summary(glm7)

anova(glm7, test = "F")

######
### PROBABILITY OF DEPARTURE
######

tagMeta2 <- read_rds('tagMeta.rds') %>%
  group_by(markerNumber) %>% 
  summarise(tagDeployStart = tagDeployStart[1]) %>%
  mutate(tagDeployStart = as.POSIXct(tagDeployStart, origin = '1970-01-01'),
         tagDeploy.jdate = as.integer(format(tagDeployStart, format = '%j')))

jdate.range <- range(flightBySunset$jdate)
depRange <- seq(jdate.range[1]-1, jdate.range[2])
numDep <- sapply(depRange, function(x){length(which(x>=flightBySunset$jdate))})

probDep <- tibble(jdate = depRange, 
                  numDep = numDep,
                  probDep = numDep/max(numDep))
probDep %>%
  ggplot(aes(jdate, probDep))+
  geom_line()

coxph.flightByWeather <- flightByWeather %>%
  filter(isDep) %>%
  mutate(year = as.factor(year(ts))) %>%
  left_join(tagMeta2, by = 'markerNumber') %>%
  filter(year(ts) == year(tagDeployStart)) %>%
  mutate(idleTime = jdate - tagDeploy.jdate,
         departure.jdate = jdate)

  test <- coxph.flightByWeather[rep(row.names(coxph.flightByWeather), coxph.flightByWeather$idleTime),] 
  test$start <- test$tagDeploy.jdate + (sequence(coxph.flightByWeather$idleTime)-1)
  test$end <- test$tagDeploy.jdate + (sequence(coxph.flightByWeather$idleTime))
  test <- test %>% mutate(event = ifelse(end == jdate, 1, 0))

coxph.flightByWeather %>%
  mutate(tagDeployStart = format(tagDeployStart, '%m-%d')) %>%
  ggplot(aes(as.Date(tagDeployStart, format = '%m-%d'), paste(age, bandsite)))+
  geom_point()+
#  geom_histogram()+
  scale_x_date()+
  facet_grid(.~year)

with(test, Surv(start, end, event))

coxph1 <- coxph(Surv(jdate) ~ age + bandsite + year + cc + wind.dir + wind.abs, data = coxph.flightByWeather)

survfit1 <- survfit(coxph1)

coxph1 <- coxph(Surv(start, end, event) ~ age + bandsite + year + cc + wind.dir + wind.abs, data = test)

summary(coxph1)


par(mfrow(c(1,1)))

plot(survfit1)
plot(as.integer(coxph1$y)[1:1193], coxph1$residuals)

length(as.integer(coxph1$residuals)[1192:1200])

with(lung, Surv(time))
heart
Surv(type = 'left')

