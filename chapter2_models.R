## Load Libraries
library(tidyverse)
library(motus)
library(stats)
library(lmtest)

# Read in flights dataframe
flightBySunset <- read_csv("flightBySunset.csv")

# Get some stats on the number of flights
flightStats <- flightBySunset %>%
  group_by(markerNumber) %>% 
  summarise(age = age[1], nflights = n()) %>% 
  group_by(age) %>% 
  summarise(a = mean(nflights), sd = sd(nflights), nflights = sum(nflights), n = n())

flightStats

####
#  MODELS
####

######
### NUMBER OF FLIGHT
######

nFlights <- flightBySunset %>% 
  group_by(markerNumber) %>% 
  summarise(age = age[1], bandsite = bandsite[1], n = n())

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


# Test: Do individuals at remote sites make more nocturnal flights and
#       do adults make fewer flights than hatch-years and show no difference between sites? 

glm1 <- nFlights %>% glm(formula = n ~ age, family = poisson)
glm2 <- nFlights %>% glm(formula = n ~ bandsite, family = poisson)

glm3 <- nFlights %>% glm(formula = n ~ age * bandsite, family = poisson)

plot(glm3)

summary(glm3)

lrtest(glm3, glm1)
lrtest(glm3, glm2)

######
### DEPARTURE TIMES (LOG)
######

dFlights <- flightBySunset %>% filter(isDep) %>% select(bySunset, age, bandsite)
  
dFlights %>%
  ggplot(aes(bandsite, log(bySunset)))+
  geom_boxplot()+
  facet_grid(age~.)

dFlights %>%
  ggplot(aes(log(bySunset))) +
  geom_histogram()+
  facet_grid(age~bandsite)

# Test:  Do individuals at remote sites depart earlier in the evening and
#        do adults depart later than hatch-years and show no difference between sites? 

glm6 <- dFlights %>% glm(formula = log(bySunset) ~ age * bandsite, family = Gamma)

plot(glm6)

summary(glm6)

anova(glm6, test = "F")

