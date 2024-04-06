####set working directory####
setwd("C:/Users/kennedy/Desktop/Cage Study R/R scripts")

####load the packages needed for the analysis#### 
install.packages("nlme")
library(nlme)

library(ggplot2)

install.packages("ggpubr")
library(ggpubr)

install.packages("dplyr")
library(dplyr)

####read in the excel file(s)#### 
#total number of foraging trips 
foraging_trips <- read.csv(file.choose("Foraging Trips.csv"), sep = ";", header = TRUE)
#average foraging duration
mean_foraging_time <- read.csv(file.choose("Average Foraging Duration.csv"), sep = ";", header = TRUE)
#total foraging duration 
sum_foraging_time <- read.csv(file.choose("Total Foraging Duration.csv"), sep = ";", header = TRUE)


####get summary stats: mean and std####
group_by(foraging_trips, treatment) %>%
  summarise(count = n(), mean = mean("life stage", na.rm = TRUE),sd = sd("trip time", na.rm = TRUE))

group_by(mean_foraging_time, treatment) %>%
  summarise(count = n(), mean = mean("life stage", na.rm = TRUE),sd = sd("trip time", na.rm = TRUE))

group_by(sum_foraging_time, treatment) %>%
  summarise(count = n(), mean = mean("life stage", na.rm = TRUE),sd = sd("trip time", na.rm = TRUE))

####visualize the data####
boxplot(Trips ~ Type, data = foraging_trips, xlab = "Life Stage",
        ylab = "Number of foraging trips", main = "Foraging Trips")

boxplot(Duration ~ Type, data = mean_foraging_time, xlab = "Life Stage",
        ylab = "Average Foraging Time", main = "Average Foraging Time")

boxplot(Duration ~ Type, data = sum_foraging_time, xlab = "Life Stage",
        ylab = "Total Foraging Time", main = "Total Foraging Time")

####Run analysis for differences between QMP+/- and Life stages: forager, nurse, newlyemerged####
##General Liner Mixed Model: Poisson Distribution##                        
lme4_model <- glm(Trips~Treatment*Type + (1|Colony), family = "poisson", data = foraging_trips)
summary(lme4_model)

lme4_model2 <- glm(Trips~Treatment+Type + (1|Colony), family = "poisson", data = foraging_trips)
summary(lme4_model2)

lme4_model3 <- glm(Duration~Treatment*Type + (1|Colony), family = "poisson", data = mean_foraging_time)
summary(lme4_model3)

lme4_model4 <- glm(Duration~Treatment+Type + (1|Colony), family = "poisson", data = mean_foraging_time)
summary(lme4_model4)

lme4_model5 <- glm(Duration~Treatment*Type + (1|Colony), family = "poisson", data = sum_foraging_time)
summary(lme4_model5)

lme4_model6 <- glm(Duration~Treatment+Type + (1|Colony), family = "poisson", data = sum_foraging_time)
summary(lme4_model6)


####Linear Mixed Effect Model: Normal Distribution####
##Foraging Trips##
nlme_model_trips <- lme(Trips~Treatment*Type, random = ~1|Colony, data = foraging_trips) 
summary(nlme_model_trips)

nlme_model_trips2 <- lme(Trips~Treatment+Type, random = ~1|Colony, data = foraging_trips) 
summary(nlme_model_trips2)

##Average Foraging Time##
nlme_model3 <- lme(Duration~Treatment*Type, random = ~1|Colony, data = mean_foraging_time) 
summary(nlme_model3)

nlme_model4 <- lme(Duration~Treatment+Type, random = ~1|Colony, data = mean_foraging_time) 
summary(nlme_model4)

##Total Foraging Time##
nlme_model5 <- lme(Duration~Treatment*Type, random = ~1|Colony, data = sum_foraging_time) 
summary(nlme_model5)

nlme_model6 <- lme(Duration~Treatment+Type, random = ~1|Colony, data = sum_foraging_time) 
summary(nlme_model6)
