## KCUR DATA ANALYSIS WITH VERSION 1 AND 2 COMBINED

library(effects)
library(ggplot2)
library(corrplot)
library(dplyr)
library(plyr)
library(tidyr)
library(survival)
library(lme4)
library(survivalAnalysis)
library(survminer)
library(coxme)
library(ggeffects)
library(frailtySurv)
library(ggfortify)
library(survAUC)
library(sjstats)

#Load cleaned files (both have already excluded participants)
KCur <- read.csv("KCur_Cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
KCur <- KCur[-1]
KCur2 <- read.csv("KCur_v2_Cleaned.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
KCur2 <- KCur2[-1]

# Bind files together
AllKCur <- rbind(KCur, KCur2)
AllKCur <- arrange(AllKCur, desc(ID), C_Trial)

# Load curiosity ratings
AllQ <- read.csv("TrivQs_Share.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AllQ <- AllQ[ -c(1, 3:4, 6:7, 9:17) ]
AllQ$C_QuestionNum <- AllQ$QuestionNum
AllQ <- AllQ[-c(1)]
AllKCur <- merge(AllKCur, AllQ, by="C_QuestionNum", all.x=TRUE)
AllKCur <- arrange(AllKCur, desc(ID), C_Trial)
AllKCur$CurZ <- scale(AllKCur$Curiosity, center=TRUE, scale=TRUE)


# Drop "Know" & "NA" trials
AllKCur <- subset(AllKCur, C_Response == "W" | C_Response == "S")
AllKCur$C_Environment <- AllKCur$C_Environment - 1 #0=UD, 1=HT
AllKCur$Observed <- ifelse(AllKCur$C_Response == "S", 1,0) #skip = 1, wait = 0 (necessary for censoring correctly!)


######### GGPLOT THEMES

theme_set(theme_bw() + theme(legend.position= c(0.7, 0.1),legend.title = element_blank(), 
                             legend.key.size = unit(1, 'cm'),
                             legend.text=element_text(size=30, margin = margin(t = 0, r = 25, b = 0, l = 0)),
                             axis.title.x = element_text(family = "Helvetica", size=30,
                                                         margin = margin(t = 20, r = 0, b = 0, l = 0)),
                             axis.title.y = element_text(family = "Helvetica", size=24, 
                                                         margin = margin(t = 0, r = 25, b = 0, l = 0)),
                             axis.text = element_text(size=30), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank()))


UD_HT_colors=c("#00BFC4", "#F8766D")

############# BEHAVIORAL ANALYSES #############

########### SURVIVAL ANALYSIS
# survival curve analysis, again excluding trials where participants waited for <= 1 second & >= 20 seconds
AllKCur_rev <- subset(AllKCur, C_WaitTime >= 1 & C_WaitTime <= 20)

# COX SURVIVAL MODEL
All_SurvCur_new2 <- coxme(Surv(C_WaitTime, Observed) ~ C_Environment * CurZ + (C_Environment |ID) + (CurZ | ID) +
                            (C_Environment |C_QuestionNum), data = AllKCur_rev)

########## CURIOSITY ANALYSIS
## WILLINGNESS TO WAIT
## How does curiosity impact people's willingness to wait in either environment? (remove trials where wait
# is less than 1, but keep all other wait times >= 20)
invlog <- function(x) { 1/(1+exp(-(x))) }

KCur4 <- subset(AllKCur, AllKCur$C_WaitTime >= 1)
KCur4$C_Response <- recode(KCur4$C_Response, "S" = 0, "W" = 1)

M2 <- glmer(C_Response ~ CurZ * C_Environment + (CurZ * C_Environment | ID) + (C_Environment|C_QuestionNum),
            data=KCur4, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
summary(M2)

pred2 <- expand.grid(C_Environment = c(0, 1), CurZ=seq(min(KCur4$CurZ), max(KCur4$CurZ), by=.1))
fit2 <- as.data.frame(predict(M2, newdata=pred2, re.form=NA))
fit2 <- cbind(pred2, fit2)
colnames(fit2)[3] <- "fit2"
fit2$Prob <- invlog(fit2$fit2)
Condse2 <- summary(M2)$coefficients[2,2]
fit2$se <- Condse2
fit2$lower <- fit2$fit2 - fit2$se*2
fit2$upper <- fit2$fit2 + fit2$se*2
fit2$lowerprob <- invlog(fit2$lower)
fit2$upperprob <- invlog(fit2$upper)
fit2$C_Environment <- as.factor(fit2$C_Environment)

ggplot(fit2, aes(x=CurZ, y=Prob, color=C_Environment, group=C_Environment)) + 
  geom_ribbon(aes(ymin=lowerprob, ymax=upperprob, fill = C_Environment), alpha=.3, color = NA) +
  geom_line(lwd=1.5) + scale_y_continuous(lim=c(0.4,1)) + labs(y = "P(Wait)", x = "Question Curiosity (Z)") +
  scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +
  scale_fill_manual(values = UD_HT_colors, labels = c("Uniform Distribution", "Heavy-Tailed Distribution"))

############### MEMORY ANALYSES
M2_New <- glmer(MC_Correct ~ CurZ * C_Environment + C_WaitTime + UD_NWaitZ + HT_NWaitZ + 
                  (CurZ + C_Environment + C_WaitTime |ID) + (C_Environment |C_QuestionNum),
                data=NewMem_df, family=binomial,
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

summary(M2_New)

##graphing the model 
M2_df <- as.data.frame(effect('CurZ*C_Environment', mod=M2_New, xlevels=list
                              (C_Environment= c(0, 1),
                                CurZ=seq(min(KCur5$CurZ), max(KCur5$CurZ), by=.1))))

Mem.Cur.fig.new <- ggplot(M2_df, aes(x=CurZ, y=fit, color=as.factor(C_Environment), group=as.factor(C_Environment))) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = as.factor(C_Environment)), alpha=.3, color=NA) +
  geom_line(lwd=3) + scale_color_manual(values=UD_HT_colors, labels=c("Uniform distribution", "Heavy-tailed distribution")) +
  scale_y_continuous(lim=c(0,1), breaks = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + labs(x="Question curiosity (z-scored)", y="P(remember)") +
  scale_fill_manual(values = UD_HT_colors, labels = c("Uniform distribution", "Heavy-tailed distribution"))

Mem.Cur.fig.new + theme(legend.background = element_rect(fill="white",
                                                         size=0.5, linetype="solid", 
                                                         colour ="black"),
                        legend.position= c(0.35, 0.25))

##Look at how information prediction error affects memory (measure of semantic surprise)
AllQ_sat <- read.csv("TrivQs_Share.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AllQ_sat <- AllQ_sat[ -c(1, 3:4, 7, 10:17) ]
AllQ_sat$C_QuestionNum <- AllQ_sat$QuestionNum
AllQ_sat <- AllQ_sat[-c(1)]
AllKCur_sat <- merge(AllKCur, AllQ_sat, by="C_QuestionNum", all.x=TRUE)
AllKCur_sat <- arrange(AllKCur_sat, desc(ID), C_Trial)
AllKCur_sat$CurZ <- scale(AllKCur_sat$Curiosity.y, center=TRUE, scale=TRUE)
AllKCur_sat$SatZ <- scale(AllKCur_sat$Satisfy, center=TRUE, scale=TRUE)
AllKCur_sat <- AllKCur_sat[ -c(20:21)]

##Compute Z-Scored info. prediction error
AllKCur_sat$PE <- AllKCur_sat$Satisfy - AllKCur_sat$Curiosity.y
AllKCur_sat$PE <- scale(AllKCur_sat$Satisfy, center=TRUE, scale=TRUE)

##Create dataframe with info. prediction errors
##Only include wait times that are bigger than one second
invlog <- function(x) { 1/(1+exp(-(x))) }
KCur5 <- subset(AllKCur_sat, AllKCur_sat$C_WaitTime >= 1)

KCur5$C_Response <- recode(KCur5$C_Response, "S" = 0, "W" = 1)

M5 <- glmer(MC_Correct ~ PE * C_Environment + (PE * C_Environment |ID) + (C_Environment|C_QuestionNum),
            data=KCur5, family=binomial,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

summary(M5)

##Plot of smoothed raw data
ggplot(KCur5) + geom_smooth(aes(x=PE, y=MC_Correct, group=as.factor(C_Environment), color=as.factor(C_Environment))) +
  scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +
  scale_fill_manual(values = UD_HT_colors, labels = c("Uniform Distribution", "Heavy-Tailed Distribution")) 

##Plot model
pred5 <- expand.grid(C_Environment = c(0, 1), PE=seq(min(KCur5$PE), max(KCur5$PE), by=.1))

fit5 <- as.data.frame(predict(M5, newdata=pred5, re.form=NA))
fit5 <- cbind(pred5, fit5)
colnames(fit5)[3] <- "fit5"
fit5$Prob <- invlog(fit5$fit5)

Condse5 <- summary(M5)$coefficients[2,2]
fit5$se <- Condse5
fit5$lower <- fit5$fit5 - fit5$se*2
fit5$upper <- fit5$fit5 + fit5$se*2
fit5$lowerprob <- invlog(fit5$lower)
fit5$upperprob <- invlog(fit5$upper)
fit5$C_Environment <- as.factor(fit5$C_Environment)

ggplot(fit5, aes(x=PE, y=Prob, color=C_Environment, group=C_Environment)) + 
  geom_ribbon(aes(ymin=lowerprob, ymax=upperprob, fill = C_Environment), alpha=.3, color = NA) +
  geom_line(lwd=1.5) + scale_y_continuous(lim=c(0,1)) + labs(y = "P(Wait)", x = "Info. Prediction Error (Z)") +
  scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +
  scale_fill_manual(values = UD_HT_colors, labels = c("Uniform Distribution", "Heavy-Tailed Distribution")) 

############# EYE-TRACKING ANALYSES #############

#Load and prepare dataframe
ShowAns_size_1s <- read.csv("~/Data/Pupil Size csv/ShowAns_size_1s.csv")
ShowAns_size_1s$environment[which(ShowAns_size_1s$environment==1)] <- 0
ShowAns_size_1s$environment[which(ShowAns_size_1s$environment==2)] <- 1

### Time course of pupil size in first second after answer is shown ###
ShowAns_size0_1s <- filter(ShowAns_size_1s, environment == 0)
c <- c(1:1000)
timeCourse_UD <- as.data.frame(replicate(1000,0))
o <- 1

for (i in c) {
  timeCourse_UD[o,1] <- mean(ShowAns_size0_1s$SizeZ[which(ShowAns_size0_1s$timepoint == i)], na.rm = TRUE)
  o <- o+1
}

o <- 1
for (i in c) {
  timeCourse_UD[o,2] <- se(ShowAns_size0_1s$SizeZ[which(ShowAns_size0_1s$timepoint == i)], na.rm = TRUE)
  o <- o+1
}

names(timeCourse_UD) <- c("SizeZ", "se")
timeCourse_UD$time <- c(1:1000)

#Plot average time course 1 second after answer appears for HT environment
ShowAns_size1_1s <- filter(ShowAns_size_1s, environment == 1)
c <- c(1:1000)
timeCourse_HT <- as.data.frame(replicate(1000,0))
o <- 1

for (i in c) {
  timeCourse_HT[o,1] <- mean(ShowAns_size1_1s$SizeZ[which(ShowAns_size1_1s$timepoint == i)], na.rm = TRUE)
  o <- o+1
}

o <- 1
for (i in c) {
  timeCourse_HT[o,2] <- se(ShowAns_size1_1s$SizeZ[which(ShowAns_size1_1s$timepoint == i)], na.rm = TRUE)
  o <- o+1
}

names(timeCourse_HT) <- c("SizeZ", "se")
timeCourse_HT$time <- c(1:1000)

#Plot both of them together
timeCourse_ShowAns.final <- ggplot() +geom_ribbon(data = timeCourse_UD, aes(x = time, ymin = SizeZ-se, ymax = SizeZ + se), fill = "#00BFC4", alpha = .3) + geom_line(data=timeCourse_UD, aes(x=time, y=SizeZ), color = "#00BFC4") + 
  geom_ribbon(data = timeCourse_HT, aes(x = time, ymin = SizeZ-se, ymax = SizeZ + se) ,fill = "#F8766D", alpha = .3) + geom_line(data=timeCourse_HT, aes(x=time, y=SizeZ), color = "#F8766D") + labs(y = "Pupil Size (Z-Scored)", x = "Time (ms)") + scale_colour_manual(values=c("#00BFC4", "#6a51a3")) +
  theme(legend.position=c(0.9, 0.9))

timeCourse_ShowAns.final <- timeCourse_ShowAns.final + geom_vline(xintercept = 0, linetype = 'dashed') 

timeCourse_ShowAns.final

### Does pupil size predict subsequent memory? ###
#Subset pupillometry data to 500 - 750 ms after the answer and run model without environment interaction
short_ShowAns <- subset(ShowAns_size_1s, timepoint >= 500 & timepoint <= 750)

ByTrial_ShowAns <- filter(short_ShowAns) %>%
  group_by(ID, trial, environment, mem, C_QuestionNum) %>%
  summarise_at(vars(SizeZ), funs(mean,sd))

#Run simple model
Mem.ByTrial_ShowAns <- glmer(mem ~ mean + (mean|ID), 
                             data = ByTrial_ShowAns,
                             family = binomial,
                             na.action = "na.omit")

summary(Mem.ByTrial_ShowAns)

#Or run model with interaction
Mem.ByTrial_ShowAns <- glmer(mem ~ mean*environment + (mean*environment|ID) + (environment|C_QuestionNum), 
                             data = ByTrial_ShowAns,
                             family = binomial,
                             na.action = "na.omit")

#Plot 
pred <- expand.grid(environment = c(0, 1), mean=seq(min(ByTrial_ShowAns$mean), max(ByTrial_ShowAns$mean), by=.1))

fit <- as.data.frame(predict(Mem.ByTrial_ShowAns, newdata=pred, re.form=NA))
fit <- cbind(pred, fit)
colnames(fit)[3] <- "fit"
fit$Prob <- invlog(fit$fit)

Condse <- summary(Mem.ByTrial_ShowAns)$coefficients[2,2]
fit$se <- Condse
fit$lower <- fit$fit - fit$se*2
fit$upper <- fit$fit + fit$se*2
fit$lowerprob <- invlog(fit$lower)
fit$upperprob <- invlog(fit$upper)

MemXPup.fig <- ggplot()  + 
  geom_ribbon(fit, mapping = aes(x=mean, y=Prob, color=as.factor(environment), group=as.factor(environment), ymin=lowerprob, ymax=upperprob, fill = as.factor(environment)), alpha=.3, color=NA) +
  geom_line(fit, mapping = aes(x=mean, y=Prob, color=as.factor(environment), group=as.factor(environment)), lwd=1.5) + scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) + scale_y_continuous(lim=c(0,1), breaks = c(0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))  + labs(x="Pupil Size (Z-Scored)", y="P(Remember)") +
  scale_fill_manual(values = UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +
  stat_summary_bin(ByTrial_ShowAns, mapping = aes(mean, mem, color = as.factor(environment), group = as.factor(environment)), fun = 'mean', bins = 100, alpha = 0.4) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.position= c(0.35, 0.2),
        legend.direction="vertical")

MemXPup.fig

### Does curiosity predict pupil size? ###
#Add curiosity rating to dataframe
AllQ$CurZ <- scale(AllQ$Curiosity, center=TRUE, scale=TRUE)
ByTrial_ShowAns <- merge(AllQ, ByTrial_ShowAns, by = "C_QuestionNum")

#Model
Cur.ByTrial_ShowAns <- lmer(mean ~ CurZ*environment + (CurZ*environment|ID), 
                            data = ByTrial_ShowAns,
                            na.action = "na.omit")

summary(Cur.ByTrial_ShowAns)

#Plot
pred <- expand.grid(environment = c(0, 1), CurZ=seq(min(ByTrial_ShowAns$CurZ), max(ByTrial_ShowAns$CurZ), by=.1))

fit <- as.data.frame(predict(Cur.ByTrial_ShowAns, newdata=pred, re.form=NA))
fit <- cbind(pred, fit)
colnames(fit)[3] <- "fit"

Condse <- summary(Cur.ByTrial_ShowAns)$coefficients[2,2]
fit$se <- Condse
fit$lower <- fit$fit - fit$se*2
fit$upper <- fit$fit + fit$se*2

PupXCur.fig <- ggplot(fit, aes(x=CurZ, y=fit, color=as.factor(environment), group=as.factor(environment))) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = as.factor(environment)), alpha=.3, color=NA) +
  geom_line(lwd=3) + scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +  
  scale_y_continuous(lim=c(-1,1)) + labs(y="Pupil Size (Z-Scored)", x="Question Curiosity") +
  stat_summary_bin(ByTrial_ShowAns, mapping = aes(mean, CurZ, color = as.factor(environment), group = as.factor(environment)), fun = 'mean', bins = 100, alpha = 0.4) +
  scale_x_continuous(lim=c(-1,3)) + 
  scale_fill_manual(values = UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) 

PupXCur.fig
####### DEMOGRPHIC INFORMATION
Demo <- read.csv("KCur_Questionnaires_ALL.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# exclude participants: 001, 002, 005, 026, 101, 108, 120, 128, 0/135, 3024
ex <- c("001", "002", "005", "026", "101", "108", "120", "128", "0", "3024", "135")
Demo <- Demo[!(Demo$ID %in% ex),]

# add in the extra participant 029 who did questionnaires in Ellen's dataset
setwd("/Users/emilylang/Dropbox/Research/ShohamyLab/KCur- Emily/Data/029")
extra <- read.csv("029_Questionnaires.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# delete everything but race, gender, age
Demo1 <- Demo[c("ID", "Age", "Sex", "Race")]
extra1 <- extra[c("ID", "Age", "Sex", "Race")]
extra1 <- extra1[-c(1:2),]
All <- rbind(Demo1, extra1)
All <- All[-c(1:2),]

All$Age <- as.numeric(All$Age)
age <- mean(All$Age)
sex <- sum(All$Sex == "Male")
indian <- sum(All$Race == "American Indian")  
africanam <- sum(All$Race == "African American")  
asian <- sum(All$Race == "Asian")  
white <- sum(All$Race == "Caucasian") 
hispanic <- sum(All$Race == "Hispanic") 
mideast <- sum(All$Race == "Middle Eastern") 
pacis <- sum(All$Race == "Pacific Islander") 
other <- sum(All$Race == "Other") 


