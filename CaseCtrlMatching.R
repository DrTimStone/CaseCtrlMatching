# Script for matching a set of cases and controls based on several co-variate factors
# Matched by different co-variates, factors are weighted
# This script will produce 192 samples, for assignment on to two epigenetic plates
# The matched covariates are:
# Age, Sex, Cigarettes (pack years), Alcohol (unit years), bmi, heartburn accumulation, PPI accumulation

library(tibble)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(gdata)
library(beepr)

HOME_DIRECTORY <- "/Users/timothystone/Desktop/CaseControlMatches/"

# Input Filelocations
# fileIn contains a comprehensive set of covariate information
# ID maps the "subject number" to the sample ids
fileIn <- paste0(HOME_DIRECTORY, "PITDataDumpMar2020.csv")
fileID <- paste0(HOME_DIRECTORY, "IDandGroups.csv")

# Output file locations
fileOut <- paste0(HOME_DIRECTORY, "CASECONTROL.xls")
filePlot <- paste0(HOME_DIRECTORY, "MATCH_PLOTS.xls")

# WEIGHTS
sexWeighting        <- 500
ageWeighting        <- 300
bmiWeighting        <- 300
weightWeighting     <- 300
heartburnWeighting  <- 300
PPIWeighting        <- 300
cigaretteWeighting  <- 100
drinkWeighting      <-  50

# Execution switches
plotSwitch      <- FALSE

# Filter switches
raceFilter      <- TRUE
salivaFilter    <- TRUE
lowGradeFilter  <- TRUE

# Constants
ACCEPTABLE_AGE_DIFFERENCE <- 11

# Initialise the looping variables
matchedControl <- c()
globalCount <- 1

# Functions

normalise <- function(x, y=100) { x * y / max(x, na.rm=T) }

# Read the ID and group data
IDandGroups <- read.csv(fileID, stringsAsFactors = F) %>% unique() %>% as_tibble()
colnames(IDandGroups) <- c("Subject.number", "FinalDiagnosis")

# Read the large sample file
samples <- read.csv(filein, stringsAsFactors = F) %>% as_tibble()
samples <- samples %>% mutate(age = ifelse(!is.na(Q_P_age) & Q_P_age > REG_ageConsent, Q_P_age, REG_ageConsent))

# Database cleanup
# Sadly database input was not screened (I didn't design it)
# People have entered text in numeric values
# Choice of unit was either kg or stone/feet, this needs to be standardised to kg / m
# The DataFrame needs to be cleaned and then coerced to a numeric format
# What follows are therefore a set of cleanup regular expressions using gsub()

samples <- samples %>% mutate(Q_P_weight_stones = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_stones = gsub("\\D+", "", Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_stones = as.numeric(Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_stones =  ifelse(Q_P_weight_stones <6 , NA, Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_pounds = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_pounds))
samples <- samples %>% mutate(Q_P_weight_pounds = gsub("\\D+", "", Q_P_weight_pounds))
samples <- samples %>% mutate(Q_P_weight_stones = as.numeric(Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_kg = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_kg))
samples <- samples %>% mutate(Q_P_weight_kg = gsub("\\D+", "", Q_P_weight_kg))
samples <- samples %>% mutate(Q_P_weight_kg = as.numeric(Q_P_weight_kg))
samples <- samples %>% mutate(Q_P_height_feet = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_feet = gsub("\\D+", "", Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_feet = as.numeric(Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_inches = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_inches))
samples <- samples %>% mutate(Q_P_height_inches = gsub("\\D+", "", Q_P_height_inches))
samples <- samples %>% mutate(Q_P_height_inches = as.numeric(Q_P_height_inches))
samples <- samples %>% mutate(Q_P_height_cm = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_cm))
samples <- samples %>% mutate(Q_P_height_cm = as.numeric(Q_P_height_cm))
samples <- samples %>% mutate(Q_P_height_cm = ifelse(Q_P_height_cm < 100, NA, Q_P_height_cm))
samples <- samples %>% mutate(height = ifelse(grepl("\\d", Q_P_height_feet) & !grepl("\\d", Q_P_height_cm),
                                              round(as.numeric(Q_P_height_feet) * 30.48 + (as.numeric(Q_P_height_inches) * 2.54),2), NA))
samples <- samples %>% mutate(height = ifelse(grepl("\\d", height), height, Q_P_height_cm))
samples <- samples %>% mutate(weight = ifelse(grepl("\\d", Q_P_weight_stones) & !grepl("\\d", Q_P_weight_kg),
                                              round(as.numeric(Q_P_weight_stones) * 6.35029 + (as.numeric(Q_P_weight_pounds) * 0.453592),2), ""))
samples <- samples %>% mutate(weight = ifelse(grepl("\\d", Q_P_weight_kg), Q_P_weight_kg, weight))
samples <- samples %>% mutate(weight = ifelse(weight == "" | weight == 0, NA, weight))
thisweight <- samples %>% select(weight) %>% unlist %>% as.numeric()
samples <- samples %>% mutate(weight = thisweight)
samples <- samples %>% mutate(weight = ifelse(weight > 1000, NA, weight))
samples <- samples %>% mutate(weight = as.numeric(weight))
samples <- samples %>% filter(!is.na(weight))
samples <- samples %>% mutate(bmi = ifelse(grepl("\\d", height) & grepl("\\d", weight), weight/(((height/100))^2), NA))
samples <- left_join(IDandGroups, samples)

# Calculation of alcohol consumption (unit years)
# The questionnaire gives general time periods and average alcohol consumption by week
# We can construct an analagous unit to the pack-year

# These are the statistifcal Bins used for taking averages from the survey data
ageBoundary <- c(17, 25, 35, 45, 59, 999)
ageRange <- c(8, 10, 10, 15)
averageUnits <- c(5, 15, 25, 35, 45)

# Extract the alcohol-related fields from the large data-file
alcoholNames <- colnames(samples[,grep("Q_ALC", colnames(samples))])[4:8]
alcoholGreps <- c("teetotal", "light drinker", "moderate drinker", "very heavy drinker")

# Initialise CumulAlc=0 in the samples DataFrame, we cannot match with NA so we need a value
samples <- samples %>% mutate(CumulAlc = 0)
ages <- select(samples, age) %>% unlist() 

# NB:  .bincode is a Base R-function that takes a numeric vector and return integer codes for the binning.
# DO NOT REMOVE THE INITAL DOT
ageBins <- .bincode (ages, ageBoundary)
unitsDrunk <- samples %>% select(CumulAlc) %>% unlist()

# Alcohol Loop
# loop through the bins, add up the accumulated drink years to get alcohol consumption to date
for (i in 1:5) {
  drinkYears <- rep(0, length(unitsDrunk))
  drinkYears[which(ageBins == i)] <- ages[which(ageBins == i)] - ageBoundary[i]
  drinkYears[which(ageBins > i)] <- ageRange[i] 
  alcoholSurvey <- select(samples, !!alcoholNames[i]) %>% unlist %>% as.character
  
  for (j in 1:5) {
    whichUnits <- grep(alcoholGreps[j], alcoholSurvey)
    unitsDrunk[whichUnits] <- unitsDrunk[whichUnits] + averageUnits[j] * drinkYears[whichUnits]
  }
}
# Add to the samples DataFrame
samples <- mutate(samples, UnitsDrunk = unitsDrunk)

# Cigarette consumption calculation, this is similar in execution to the alcohol consumption

ageBoundary <- c(15, 25, 35, 45, 59, 999)
ageRange <- c(10, 10, 10, 15)
averageCigs <- c(5, 15, 30, 50)
smokeNames <- colnames(samples[,grep("Q_SMK", colnames(samples))])[5:9]
smokeGreps <- c("10", "20", "21", "More than")
samples <- samples %>% mutate(CumulSmoke = 0)
ages <- select(samples, age) %>% unlist() 
ageBins <- .bincode (ages, ageBoundary) # .bincode AGAIN, KEEP THE DOT
cigsSmoked <- samples %>% select(CumulSmoke) %>% unlist()

# Loop through the bins, as per alcohol
for (i in 1:5) {
  smokeYears <- rep(0, length(cigsSmoked))
  smokeYears[which(ageBins == i)] <- ages[which(ageBins == i)] - ageBoundary[i]
  smokeYears[which(ageBins > i)] <- ageRange[i] 
  smokeSurvey <- select(samples, !!smokeNames[i]) %>% unlist %>% as.character
  
  for (j in 1:4) {
    whichSmoke <- grep(smokeGreps[j], smokeSurvey)
    cigsSmoked[whichSmoke] <- cigsSmoked[whichSmoke] + averageCigs[j] * smokeYears[whichSmoke]
  }
}

samples <- mutate(samples, CigsSmoked = cigsSmoked)

# PPI Consumption

samples <- mutate(samples, PPIfreq = as.character(Q_HB_ppi_medication_frequency))
samples <- mutate(samples, PPIfreq = gsub("(.+) $", "\\1", PPIfreq, perl=T))

# This section turns PPI estimates into numbers for yearly consumption, hence daily = 365
samples <- mutate(samples, PPIfreq = case_when(
  PPIfreq == "Daily" ~ 365,
  PPIfreq == "Few times a month" ~ 24,
  PPIfreq == "Few times a week" ~ 150,
  PPIfreq == "Few times a year" ~ 6,
  PPIfreq == "Never" ~ 0,
  PPIfreq == "" ~ 0))

# This takes the time that people began taking PPIs
samples <- mutate(samples, PPIbeg = as.character(Q_HB_ppi_medication_begin))
samples <- mutate(samples, PPIbeg = case_when(
  PPIbeg == "1 to 5 years" ~ 2.5,
  PPIbeg == "10 to 20 years" ~ 15,
  PPIbeg == "5 to 10 years" ~ 7.5,
  PPIbeg == "6 months to 1 years" ~ 0.75,
  PPIbeg == "Less than 6 months" ~ 0.25,
  PPIbeg == "More than 20 years" ~ 30,
  PPIbeg == "" ~ 0))
samples <- mutate(samples, PPIbeg = ifelse(is.na(PPIbeg), 0, PPIbeg))

# If we don't have PPI beginning frequency
# Estimate that if someone previous took PPIs they did it for a fifth of their adult life
samples <- mutate(samples, adultLife = age - 18)
samples <- mutate(samples, adultDoseLife = adultLife / 5)
samples <- mutate(samples, PPIprev = as.character(Q_HB_previous_ppi_medication_frequency))

# This section processess people that USED to take PPI
samples <- mutate(samples, PPIprev = case_when(
  PPIprev == "Daily" ~ 365,
  PPIprev == "Few times a month" ~ 24,
  PPIprev == "Few times a week" ~ 150,
  PPIprev == "Few times a year" ~ 6,
  PPIprev == "Never" ~ 0,
  PPIprev == "" ~ 0))
samples <- mutate(samples, ifelse(is.na(PPIprev), 0, PPIprev))
samples <- mutate(samples, PPIprev = ifelse(PPIfreq == 0, PPIprev, 0))

# Final PPI estimation value, a new PPI column is created for "PPI years"
samples <- mutate(samples, PPI = (PPIfreq * PPIbeg) + (PPIprev * adultDoseLife))

# Heartburn calculation
# Heartburn does not have a beginning field, so we use only current frequency

samples <- mutate(samples, hb = as.character(Q_HB_heartburn_previous_frequency))
samples <- mutate(samples, hb = case_when(
  hb == "Daily" ~ 365,
  hb == "Few times a month" ~ 24,
  hb == "Few times a week" ~ 150,
  hb == "Few times a year" ~ 6,
  hb == "Never" ~ 0,
  hb == "" ~ 0))

diagnosisRef <- grep("PAT_highest_grade_diagnosis_ever$", colnames(samples), perl=T)
diagnosis <- colnames(samples)[diagnosisRef]
diagnosisNo <- diagnosisRef

groups <- samples[,diagnosisRef] %>% unlist()

# The different degrees of cancer will be placed in a "case" group
# The different controls in a "control" group
# The CaseControlVector will be inserted into samples

# Initialise Group Vector with NA values
caseControlVector <- rep(NA, length(groups))
caseControlLevels <- samples[,diagnosisRef] %>% unlist %>% levels

case <-     caseControlLevels[3:6]
control <-  caseControlLevels[7:8]
unknown <-  caseControlLevels[1:2]
other <-    caseControlLevels[9]

# Turn the caseControlVector into a character vector (and ultimately into.a factor)
caseControlVector[groups %in% case]     <-  "Case"
caseControlVector[groups %in% control]  <-  "Control"
caseControlVector[groups %in% unknown]  <-  "DontKnow"
caseControlVector[groups %in% other]    <-  "Other"

# Insert the group vector into samples
samples <- mutate(samples, Group = caseControlVector)

# everyone = everyone eligible for case-control matching
everyone  <- filter(samples, !is.na(FinalDiagnosis))

# Remove samples with missing gender OR age information
everyone <-  filter(everyone, is.na(Gender) | is.na(age))

# Apply Filters

# Race filter
if (raceFilter) {
  #Race Filter
  everyone <- filter(everyone, grepl("White", PINF_ethnicity) | grepl("White", Q_P_ethnicity) )
}

# Low-Grade Dysplasia filter
if (lowGradeFilter) {
  #Low-grade Dysplasia Filter
  everyone <- filter(everyone, diagnosis != "Low Dys")
}

# Saliva Sample Filter
if (salivaFilter) {
  # Only people with saliva samples
  everyone <- filter(everyone, grepl("aliva", QS_samples_collected))
}


# Shorten the diagnosis names (makes for easier output)
everyone  <- mutate(everyone, diagnosis = as.character(PAT_highest_grade_diagnosis_ever))
everyone  <- mutate(everyone, diagnosis = as.character(diagnosis))
everyone  <- mutate(everyone, diagnosis = 
                      case_when(
                        diagnosis == "High Grade Dysplastic Barrett's oesophagus" ~ "High Dys",
                        diagnosis == "Low Grade Dysplastic Barrett's oesophagus" ~ "Low Dys",
                        diagnosis == "Intramucosal Oesophageal Adenocarcinoma" ~ "Intramuc. OAC",
                        diagnosis == "Invasive Oesophageal Adenocarcinoma" ~ "Invasive OAC",
                        diagnosis == "Non dysplastic Barrett's oesophagus" ~ "NDBE",
                        grepl("now", diagnosis) ~ "Don't know",
                        diagnosis== "Normal (no Barrett's)" ~ "Normal",
                        grepl("Other ", diagnosis) ~ "Other")) 

# Turn the diagnoses from a character vector into a factor, assign levels 
everyone <- mutate(everyone, diagnosis = factor(diagnosis, 
                                                levels = c("Normal", "NDBE", "Low Dys", "High Dys",
                                                           "Intramuc. OAC", "Invasive OAC", "DontKnow", "Other")))

cases     <- filter(everyone, Group == "Case")
controls  <- filter(everyone, Group == "Controls")

# Optional plotting of covariate data
if (plotSwitch) {
  jpeg(filePlot)
  agevector <- c()
  ggplot(data = everyone, aes(x=diagnosis, y = CigsSmoked, col=diagnosis)) +  geom_boxplot(notch = F)
  
  ggplot(data = everyone, aes(x=diagnosis, y = UnitsDrunk, col=diagnosis)) +  geom_boxplot(notch = F)
  
  ggplot(data = agevector, aes(x=Age, y=No.of.Younger.People, col=Sex)) + geom_line(size=2)
  
  ggplot(data = everyone, aes(x=diagnosis, y = weight, col=diagnosis)) + geom_violin()
  dev.off()
  
}

## Normalise all the quantities, create new DataFrame quantities of the normalised quantities

everyone <- everyone %>% mutate(Nweight = normalise (weight, weightWeighting))
everyone <- everyone %>% mutate(NcigsSmoked = normalise (cigsSmoked, cigaretteWeighting))
everyone <- everyone %>% mutate(Nbmi = normalise (bmi, bmiWeighting))
everyone <- everyone %>% mutate(NunitsDrunk = normalise (unitsDrunk, drinkWeighting))
everyone <- everyone %>% mutate(Nage = normalise (age, ageWeighting))
everyone <- everyone %>% mutate(NPPI = normalise(PPI, PPIWeighting))
everyone <- everyone %>% mutate(Nhb = normalise(hb, heartburnWeighting))
everyone <- mutate(everyone, nGenderVec = normalise(as.numeric(factor(Gender))-1, sexWeighting))
everyone <- mutate(everyone, nIMvector = ifelse(FinalDiagnosis == "IM" | FinalDiagnosis == "Cancer", imcWeighting, 0))
everyone <- mutate(everyone, clusterinfo = paste(Subject.number,round(age,0), paste(round(weight,0), "kg", sep=""),
                                                 CigsSmoked,FinalDiagnosis,sep="_"))
everyone <- mutate(everyone, weight = round(weight, 1))
everyone <- mutate(everyone, bmi = round(bmi, 1))

everyoneCase <- filter (FinalDiagnosis %in% c("Cancer", "HG")) %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, 
                                                                          Nage, Nweight, Nbmi, NcigsSmoked, NPPI, Nhb, NunitsDrunk) 

everyoneControls <- filter (FinalDiagnosis %in% c("Normal", "HV", "IM")) %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, 
                                                                                    Nage, Nweight, Nbmi, NcigsSmoked, NPPI, Nhb, NunitsDrunk)

normalizedQuantitiesRef <- grep("N\\w+",colnames(everyoneCase))

# A status column is added in order to flag unacceptable age matches or gender mismatches
everyoneCase <- sample_frac(everyoneCase)
everyoneCase <- mutate(everyoneCase, Status = "")
everyoneControls <- sample_frac(everyoneControls)
everyoneControls <- mutate(everyoneControls, Status = "")

# initialise output dataframe
everyMatch <- c()

# THE DATA HAS BEEN PREPARED
# THIS IS THE CASE CONTROL MATCH 

# THEIS IS THE MAIN CASE-CONTROL MATCHING LOOP
# We know in advanced that the controls outnumber the cases, so we loop through the cases and stop
# The DataFrame everyone is reconstructed from the Cases and Control sub-dataframes on each loop
# But gets smaller on each iteration

for (i in 1:nrow(neveryoneCase)) {
  
  everyone <- rbind(everyoneCase, everyoneControls)
  
  # The dist function creates an all against all match
  distancePairs <- dist(everyoneCopy[, normalizedQuantitiesRef]) %>% as.matrix()
  
  # We need to ignore cancer-cancer matches, and self-self matches
  # This is achieved by giving them undesirable values, i.e. the max value * 1000
  distancePairsMax <- max(distancePairs)[1][[1]]
  
  # Stop the subjects matching to themselves, make the central diagonal huge
  distancePairs[row(distancePairs) == col(distancePairs)]  <- distancePairsMax * 1000
  
  # Stop any Cancer-cancer or Control-control matching, make those quadrants huge
  distancePairs[col(distancePairs) <= ncase] <- distancePairsMax * 1000
  distancePairs[row(distancePairs) > ncase] <- distancePairsMax * 1000
  
  # the matched best pair now is the shortest distance, note that there are two best matches, due to matrix symmetry
  # Select the first one, hence the indexing on the which() function
  matchedIndices <- which(distancePairs == min(distancePairs), arr.ind = T)[1,]
  
  # These are the co-ordinate references of the best match
  bestCase    <- matchedIndices[1]
  bestControl <- matchedIndices[2]
  
  # Now check the matches for an unacceptable age match or a gender mismatch and flag it
  if (abs(everyoneCase[bestCase,"age"] - everyoneControls[bestControl,"age"]) > ACCEPTABLE_AGE_DIFFERENCE) {
    everyoneCase[bestCase, "status"] <- "AGE MISMATCH"
    everyoneControl[bestControl, "status"] <- "AGE MISMATCH"
  }
  if (everyoneCase[bestCase,"Gender"] != everyoneControls[bestControl,"Gender"]) {
    everyoneCase[bestCase, "status"] <- paste(everyoneCase[bestCase, "status"], "GENDER MISMATCH")
    everyoneCase[bestControl, "status"] <- paste(everyoneCase[bestControl, "status"], "GENDER MISMATCH")
  }
  
  # Augment the results vector wiht the best matched pair
  everyMatch <- rbind(everyMatch, everyoneCase[bestCase, ], everyoneControl[bestControl,])
  
  # Take out the matches from the Case and Control Matrices
  everyoneCase <- everyoneCase[-bestCase,]
  everyoneControls <- everyoneControl[-bestControl,]
  
}

# Write the output 
write.csv(everyMatch, fileOut, quote=F, row.names=F)

cat("The Cases have everyone been Matched")
beep(8)