##  Divergence analysis from EyeTrackingR
#  PG - Oct, 2018
#
#  THIS VERSION is for MVOR 3 - 1400 msec ISI with BOTH target and comp in array
# Currently being updatedd for 2250 msec ISI -KS

library(reshape2)
library(gtools) 
library(lme4)
library(plyr)
library(lattice)
library(car)
library(HH)
library(agricolae) 
library(multcomp) 
library(MASS)
library(heplots) 
library (coin) 
library(ggplot2)
library(Hmisc)
library(lmerTest)
library(LMERConvenienceFunctions)
library(languageR)
library(doBy)
library(RColorBrewer)
library("eyetrackingR")
library(pbapply)
library(ggsn)

#source ("eyeRTrackLossCreateColumn.r")
source("MVOR-Scripts-fun3_22_19.r") # PG-defined functions for script

#MV2A <- read.csv(file = "MV2A_EyetrackingR_Kids.csv", header=TRUE)
#MV2A_SMP <- read.csv(file = "MV2_SMP_EXP_ONLY_1400_TSTAMP.csv", header=TRUE)
#MV2A_SMP <- read.csv(file = "MV2_SMP_FILLER_1400_TSTAMP.csv", header=TRUE)

MV2_SMP_FILLER <- read.csv(file = "MV2_SMP_2250ms_FILLER_TSTAMP_ADULTS.csv", header = T)
MV2_SMP_EXP <- read.csv(file = "MVOR2_2250ms_Adults_TSTAMP.csv", header = T)

##  Subjects determined by prior analysis to have too many trials above 50% loss 
##  when Non-AOI looks are excluded
##  Still need to run trackloss to exclude additional 49 trials above 50%
##  No other subject has more than 6 trials above 50%
##  This is for ADULTS in the 1400 msec ISI condition!
#InReadData <- subset(MV2A_SMP, Subject != xxx)
#InReadData <- subset(InReadData, Subject != xxx)

###Creating a column that specifies whether the eye tracker lost the eye for a given sample
MV2_SMP_EXP$SampleLost <- 999

#Attempt at lapply use in this setup
#InReadData$SampleLost <- lapply( InReadData, function (x) { (InReadData$L.Raw.X..px. == 0) & (InReadData$L.POR.X..px.==0)} )

#REPLACE w. dplyr call
for (i in 1:nrow(MV2_SMP_EXP)) {
  #Target
  if ((MV2_SMP_EXP[i, c("L.POR.X..px.")] == 0) & (MV2_SMP_EXP[i, c("L.POR.Y..px.")] == 0)) {
    MV2_SMP_EXP[i, c("SampleLost")] <- 1
  } else {
    MV2_SMP_EXP[i, c("SampleLost")] <- 0
  }
  if ((i %% 1000) == 0) {
    print(i)
  }
}

write.csv(MV2_SMP_EXP, file = "MV2A_Adults_2250wSampleLostColumn.csv")
MV2_SMP_EXP_SL <- read.csv(file = "MV2A_Adults_2250wSampleLostColumn.csv", header=TRUE)

unique(MV2_SMP_EXP_SL$Subject)

####create file for eyetrackingR
MV2A_ETR_ADULT <- make_eyetrackingr_data(MV2_SMP_EXP, 
                                         participant_column = "Subject",
                                         trial_column = "Trial",
                                         time_column = "SampleTime",
                                         trackloss_column = "SampleLost",
                                         aoi_columns = c('Target','Complement', 'DiffExemplar', 'Distractor'),
                                         treat_non_aoi_looks_as_missing = TRUE)

###########################################################################################
#######Look at track loss, and remove participants and trials with too much track loss#####
###########################################################################################
#Create a file removing any time points after 2000ms 
MV2A_ETR_2000 <- subset_by_window(MV2A_ETR_ADULT,
                                  window_start_time = 0, 
                                  window_end_time = 2000, 
                                  rezero = TRUE, 
                                  remove = TRUE)

#How many subjects
unique(MV2A_ETR_2000[, c('Subject')])

###Look at track loss data (data where eyetracker lost the eye)
trackloss <- trackloss_analysis(data = MV2A_ETR_2000)
trackloss
write.csv(trackloss, file = "MV2A_Adults_2250_TracklossNoNonAOI.csv")

#trackloss <- read.csv("MV2A_Adults_TracklossNoNonAOI.csv")

#Graph trackloss by subject
meanTLoss_by_Subj <- trackloss %>% 
  group_by(Subject) %>% 
  summarize(averaged.TLoss = mean(TracklossForParticipant))

#Graph samples remaining/present per subject
MeanSamplesPresent  <- trackloss %>% 
  group_by(Subject) %>% 
  summarize(averaged.Samples = mean(Samples))

describe(MeanSamplesPresent)
describe(meanTLoss_by_Subj)

###Remove subjects and trials with over a proportion certain proportion of track loss
MV2A_Clean <- clean_by_trackloss(MV2A_ETR_2000, participant_prop_thresh = 0.5, 
                                          trial_prop_thresh = .50)

unique(MV2A_Clean$Subject)
########################################################
###     Descriptives and figure time!            #######
########################################################
#
# Create proportion of looks across two AOIs for each of the features, parts and 
#   target-only ("fillers") conditions, in time bins of 25 msec
#
TimeBins <- make_time_sequence_data(MV2A_Clean, 
                                    time_bin_size = 25, 
                                    predictor_columns = c("Condition", "TrialType"),
                                    aois = c("Target", "Complement", "DiffExemplar", "Distractor"),
                                    summarize_by = "Subject"
                                    )
colnames(TimeBins)
TimeBins$AOI<-as.factor(TimeBins$AOI)

####Remove NaN's from proportion
TB_NoNA <- TimeBins[TimeBins$Prop != 'NaN', ]
#TB_NoNA

# How many NAs?
TB_NAs <- TimeBins[TimeBins$Prop == 'NaN', ]
#TB_NAs

###Create df with only Experimental Trials
TB_Exp <- TB_NoNA[TB_NoNA$TrialType == 'Experimental', ]

# Within that df, restric to Parts only
TB_Exp_Parts <- TB_Exp[TB_Exp$Condition == 'Parts', ]

# Within that df, restric to Parts only
TB_Exp_Feat <- TB_Exp[TB_Exp$Condition == 'Features', ]

###Create file with only Filler Trials
# And divide into parts & features
TB_Filler <- TB_NoNA[TB_NoNA$TrialType == 'Filler', ]
TB_Filler_Parts <- TB_Filler[TB_Filler$Condition == "Parts", ]
TB_Filler_Feat <- TB_Filler[TB_Filler$Condition == "Features", ]

######################
##  PARTS:
#####################

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_TvsC <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'Complement') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_TvsC)

############################################
######## Target vs, DiffExempar analysis####
############################################

###Create file with only target and Diff Exemp
LookProp_Exp_Parts_TvsDE <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'DiffExemplar') , ]

myClusterTTest_Compare (LookProp_Exp_Parts_TvsDE)
#myCluster_LMER_Compare(LookProp_Exp_Parts_TvsDE)

############################################
######## Target vs. Distractor analysis#####
############################################
LookProp_Exp_Parts_TvsDist <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_TvsDist)
#myClusterCompare (LookProp_Exp_Parts_TvsDist)

############################################
######## Complement vs, DiffExempar analysis####
############################################
LookProp_Exp_Parts_CompvsDE <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Complement' | TB_Exp_Parts$AOI == 'DiffExemplar') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_CompvsDE)
#myCluster_LMER_Compare(LookProp_Exp_Parts_CompvsDE)

############################################
######## Complement vs. Distractor analysis####
############################################
LookProp_Exp_Parts_CompvsDist <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Complement' | TB_Exp_Parts$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_CompvsDist)
#myCluster_LMER_Compare(LookProp_Exp_Parts_CompvsDE)


#########################
###  FEATURES
########################

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsC <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'Complement') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsC)

#LookProp_Exp_Feat_TvsC$TargetRef<- relevel(LookProp_Exp_Feat_TvsC$AOI, ref = "Target")

############################################
######## Target vs, DiffExempar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsDE <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'DiffExemplar') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsDE)

############################################
######## Target vs. Distractor analysis#####
############################################

LookProp_Exp_Feat_TvsDist <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsDist)

############################################
######## Complement vs. Diff Exemplar  #####
############################################

LookProp_Exp_CompvsDiffExem <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Complement' | TB_Exp_Feat$AOI == 'DiffExemplar') , ]
myClusterTTest_Compare (LookProp_Exp_CompvsDiffExem)
#####################
#FixProp <- ddply(TB_Exp, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp
nrow(FixProp)

############################################
######## Complement vs. Distractor  #####
############################################

LookProp_Feat_CompvsDist <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Complement' | TB_Exp_Feat$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Feat_CompvsDist)
#####################
#FixProp <- ddply(TB_Exp, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp
nrow(FixProp)
###################################
###Parts x AOI Figure##############
###################################

source("MVOR-Scripts-fun3_22_19.r") # PG-defined functions for script
FixProp <- ddply(TB_Exp, c("Condition", "AOI", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Parts <- FixProp[FixProp$Condition == 'Parts', ]
rownames(FixProp_Parts) <- NULL

###Reorder bars### 
FixProp_Parts <- FixProp_Parts[c(241:320, 1:80, 81:160, 161:240), ]
FixProp_Parts$AOI <- factor(FixProp_Parts$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))
myGGplot (FixProp_Parts, "Parts Prime")

###################################
###  Features x AOI Figure  #######
###################################

source("MVOR-Scripts-fun3_22_19.r") # PG-defined functions for script
FixProp <- ddply(TB_Exp, c("Condition", "AOI", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Feat <- FixProp[FixProp$Condition == 'Features', ]
rownames(FixProp_Feat) <- NULL

###Reorder bars###
FixProp_Feat <- FixProp_Feat[c(241:320, 1:80, 81:160, 161:240), ]
FixProp_Feat$AOI <- factor(FixProp_Feat$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))
myGGplot(FixProp_Feat, "Feature Prime")

#######################################################
##  FILLER (TARGET ONLY) TRIALS:
## This is for Parts and Features! Separate trials
######################################################
source("MVOR-Scripts-fun.r") # PG-defined functions for script

Test_Filler_Parts <- TB_Filler_Parts[(TB_Filler_Parts$AOI == 'Target' | TB_Filler_Parts$AOI == 'DiffExemplar') , ]
Test_Filler_Feat <- TB_Filler_Feat[(TB_Filler_Feat$AOI == 'Target' | TB_Filler_Feat$AOI == 'DiffExemplar') , ]

# This step averages across the 3 distractors so that the T-test can look at only one distractor level
# Does not quite work!  PG 12/14/18 - need to actually transform to average this...
TB_Filler$AOI <- ifelse((TB_Filler$AOI == 'Target'), 'Target', 'Distractor')
# This should combine the 3 distractors into one but it does not
TB_Filler_Ave <- ddply(TB_Filler, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))

myClusterTTest_Compare (Test_Filler_Parts)
myClusterTTest_Compare (Test_Filler_Feat)

#Re-code to compare Parts to Features across trials (Target-Only trials)
TB_Target_Only <- TB_Filler[TB_Filler$AOI == "Target", ]
TB_Filler$AOI_Comp <- ifelse((TB_Filler$Condition == 'Parts'), 'Target', 'Distractor')

###Make Variable so I can average across all the distractors
TB_Filler$AOI_Filler <- ifelse((TB_Filler$AOI == 'Target'), 'Target', 'Distractor')

###Fixation proportion###
#FixProp_Filler

############################################
# Compare Parts and Features on fix prop
############################################
Test_Filler_AllTarg <- TB_Filler[(TB_Filler$AOI == 'Target'), ]
myClusterTTest_CompareBtwnConditions (Test_Filler_AllTarg)

###Reorder bars###
FixProp_Filler <- ddply(TB_Filler, c("Condition", "AOI_Filler", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Filler <- FixProp_Filler[c(193:256, 65:128, 129:192, 1:64), ]

FixProp_Filler$ArrayObject[1:64] <- "PartsTarget"
FixProp_Filler$ArrayObject[65:128] <- "FeatTarget"
FixProp_Filler$ArrayObject[129:192] <- "PartsDistr"
FixProp_Filler$ArrayObject[193:256] <- "FeatDistr"

source("MVOR-Scripts-fun.r") # PG-defined functions for script
FixProp_Filler$AOI_Filler <- factor(FixProp_Filler$ArrayObject, levels = c("PartsTarget", "FeatTarget", "PartsDistr", "FeatDistr"))
myFillerGGplot(FixProp_Filler, "Target Only")