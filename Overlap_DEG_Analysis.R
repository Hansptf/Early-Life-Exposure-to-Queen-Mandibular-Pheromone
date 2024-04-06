#set working directory for all DEG files for QMP/Cage Study
setwd <- ("C:/Users/apk6171/OneDrive - The Pennsylvania State University/Cage Study/")

#load necessary packages for analysis 
install.packages("dplyr")
library(dplyr)

##### Antennae #####
#all files for Antennae (AT)
#control vs. QMP+ (cage effect)
cage_1_day <- read.csv("New_Foragers 24day AT QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_7_day <- read.csv("Nur_Foragers 24day AT QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_19_day <- read.csv("Foragers-Foragers 24 day AT QMP+ vs Control Forager-DEGs.csv", header = TRUE, sep = ",")

#overlap between DEG lists
cage_AT <- Reduce(intersect,list(cage_1_day$Geneid, cage_7_day$Geneid, cage_19_day$Geneid))  
#save the output as the shared geneIDs for cage effect for AT  
write.csv(cage_AT, "cage all shared DEGs AT.csv", row.names = TRUE)
head(cage_AT)
summary(cage_AT)

#overlap between 1 vs 7 day
cage_AT1 <- Reduce(intersect, list(cage_1_day$Geneid, cage_7_day$Geneid))
write.csv(cage_AT1, "cage 1 vs 7 day shared DEGs AT.csv", row.names = TRUE)
head(cage_AT1)
summary(cage_AT1)

#overlap between 1 vs 19 day
cage_AT2 <- Reduce(intersect, list(cage_1_day$Geneid, cage_19_day$Geneid))
write.csv(cage_AT2, "cage 1 vs 19 day shared DEGs AT.csv", row.names = TRUE)
head(cage_AT2)
summary(cage_AT2)


#QMP- vs. QMP+ (treatment effect)  
treatment_1_day <- read.csv("New_Foragers 24day AT QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
treatment_7_day <- read.csv ("Nur_Foragers 24day AT QMP- vs QMP+ DEGs.csv", header = TRUE, sep =",")
  
treatment_19_day <- read.csv("Foragers-Foragers 24 day AT QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")

#overlap between DEG lists
treatment_AT <- Reduce(intersect,list(treatment_1_day$Geneid, treatment_7_day$Geneid, treatment_19_day$Geneid))  
#save the output as the shared geneIDs for QMP treatment for AT
write.csv(treatment_AT, "treatment all shared DEGs AT.csv ",row.names = TRUE) 
head(treatment_AT)
summary(treatment_AT)

treatment_AT1 <- Reuce(intersect, list(treatment_1_day$Geneid, treatment_19_day$Geneid))
write.csv(treatment_AT1, "treatment 1 vs 19 shared DEGs AT.csv ",row.names = TRUE) 
head(treatment_AT1)
summary(treatment_AT1)

treatment_AT2 <- Reduce(intersect, list(treatment_1_day$Geneid, treatment_7_day$Geneid))
write.csv(treatment_AT2, "treatment 1 vs 7 shared DEGs AT.csv ",row.names = TRUE) 
head(treatment_AT2)
summary(treatment_AT2)


##### Antennal Lobes #####  
#all files for Antennal Lobes (AL)
#control vs. QMP+ (cage effect)
cage_1_day <- read.csv("New_Foragers 24day AL QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_7_day <- read.csv("Nur_Foragers 24day AL QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_19_day <- read.csv("Foragers-Foragers 24 day AL QMP+ vs Control Forage-DEGs.csv", header = TRUE, sep = ",")

#overlap between DEG lists
cage_AL <- Reduce(intersect,list(cage_1_day$Geneid, cage_7_day$Geneid, cage_19_day$Geneid))  
#save the output as the shared geneIDs for cage effect for AL   
write.csv(cage_AL, "cage all shared DEGs AL.csv", row.names = TRUE) 
head(cage_AL)
summary(cage_AL)

cage_AL1 <- Reduce(intersect, list(cage_1_day$Geneid, cage_19_day$Geneid))
write.csv(cage_AL1, "cage 1 vs 19 shared DEGs AL.csv", row.names = TRUE) 
head(cage_AL1)
summary(cage_AL1)

cage_AL2 <- Reduce(intersect, list(cage_1_day$Geneid, cage_7_day$Geneid))
write.csv(cage_AL2, "cage 1 vs 7 shared DEGs AL.csv", row.names = TRUE) 
head(cage_AL2)
summary(cage_AL2)

#QMP- vs. QMP+ (treatment effect)  
treatment_1_day <- read.csv("New_Foragers 24day AL QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
treatment_7_day <- read.csv("Nur_Foragers 24day AL QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
treatment_19_day <- read.csv("Foragers-Foragers 24 day AL QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")

#overlap between DEG lists
treatment_AL <- Reduce(intersect,list(treatment_1_day$Geneid, treatment_7_day$Geneid, treatment_19_day$Geneid))  
#save the output as the shared geneIDs for QMP treatment for AL
write.csv(treatment_AL, "treatment all shared DEGs AL.csv", row.names = TRUE)
head(treatment_AL)
summary(treatment_AL)

treatment_AL1 <- Reduce(intersect, list(treatment_1_day$Geneid, treatment_19_day$Geneid))
write.csv(treatment_AL1, "treatment 1 vs 19 shared DEGs AL.csv", row.names = TRUE)
head(treatment_AL1)
summary(treatment_AL1)

treatment_AL2 <- Reduce(intersect, list(treatment_1_day$Geneid, treatment_7_day$Geneid))
write.csv(treatment_AL2, "treatment 1 vs 7 shared DEGs AL.csv", row.names = TRUE)
head(treatment_AL2)
summary(treatment_AL2)


##### Mushroom Bodies #####  
#all files for the Mushroom bodies (MB)
#control vs. QMP+ (cage effect)
cage_1_day <- read.csv("New_Foragers 24day MB QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_7_day <- read.csv("Nur_Foragers 24day MB QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
  
cage_19_day <- read.csv("Foragers-Foragers 24 day MB QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")

#overlap between DEG lists
cage_MB <- Reduce(intersect,list(cage_1_day$Geneid, cage_7_day$Geneid, cage_19_day$Geneid))
#save the output of shared geneIDs as overlap between cage effects in MB
write.csv(cage_MB, "cage all shared DEGs MB.csv", row.names = TRUE)   
head(cage_MB)
summary(cage_MB)

cage_MB1 <- Reduce(intersect, list(cage_1_day$Geneid, cage_19_day$Geneid))
write.csv(cage_MB1, "cage 1 vs 19 shared DEGs MB.csv", row.names = TRUE)
head(cage_MB1)
summary(cage_MB1)

#QMP- vs. QMP+ (treatment effect)  
treatment_1_day <- read.csv("New_Foragers 24day MB QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
treatment_7_day <- read.csv("Nur_Foragers 24day MB QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
treatment_19_day <- read.csv("Foragers-Foragers 24 day MB QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
  
#overlap between DEG lists
treatment_MB <- Reduce(intersect,list(treatment_1_day$Geneid, treatment_7_day$Geneid, treatment_19_day$Geneid))
#save the output of shared geneIDs as overlap between QMP treatments in MB
write.csv(treatment_MB, "treatment all shared DEGs MB.csv", row.names = TRUE)
head(treatment_MB)
summary(treatment_MB)

treatment_MB1 <- Reduce(intersect, list(treatment_1_day$Geneid, treatment_19_day$Geneid))
write.csv(treatment_MB1, "treatment 1 vs 19 shared DEGs MB.csv", row.names = TRUE)
head(treatment_MB1)
summary(treatment_MB1)


###### Testing the overlap between treatment groups #####

### Newly emerged ###
## treatment: QMP+ vs QMP- ##
treatment_1day_AT <- read.csv("New_Foragers 24day AT QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
treatment_1day_AL <- read.csv("New_Foragers 24day AL QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
treatment_1day_MB <- read.csv("New_Foragers 24day MB QMP- vs QMP+ DEGs.csv", header = TRUE, sep = ",")
# check the geneID overlap between each DEG list for each tissue #
treatment_1day <- Reduce(intersect, list(treatment_1day_AT$Geneid, treatment_1day_AL$Geneid, treatment_1day_MB$Geneid))
write.csv(treatment_1day, "treatment all shared DEGs 1 day.csv", row.names = TRUE)
head(treatment_1day)
summary(treatment_1day)
## cage: control foragers vs QMP+ ##
cage_1day_AT <- read.csv("New_Foragers 24day AT QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
cage_1day_AL <- read.csv("New_Foragers 24day AL QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
cage_1day_MB <- read.csv("New_Foragers 24day MB QMP+ vs Control Foragers-DEGs.csv", header = TRUE, sep = ",")
cage_1day <- Reduce(intersect, list(cage_1day_AT$Geneid, cage_1day_AL$Geneid, cage_1day_MB$Geneid))
write.csv(cage_1day, "cage all shared DEGs 1 day.csv", row.names = TRUE)
head(cage_1day)
summary(cage_1day)
