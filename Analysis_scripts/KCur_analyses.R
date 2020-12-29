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

######## SURVIVAL ANALYSIS
# survival curve analysis, again excluding trials where participants waited for <= 1 second & >= 20 seconds
AllKCur_rev <- subset(AllKCur, C_WaitTime >= 1 & C_WaitTime <= 20)

# KAPLAN MEIER FORM: no covariates/regressors/random effects; just with environment as a predictor
M1 <- survfit(Surv(C_WaitTime, Observed) ~ C_Environment, data = AllKCur_rev)
ggsurvplot(M1, data = AllKCur_rev, conf.int = TRUE, censor = FALSE, size = 1, palette = c("#00BFC4", "#F8766D"), 
           legend = "right", legend.labs = c("UD", "HT"), xlab="Wait Time (s)",
           ggtheme = theme_bw() + theme(legend.title = element_blank(), 
                                        legend.key.size = unit(1.5, 'lines'),
                                        legend.text=element_text(size=18, margin = margin(t = 0, r = 25, b = 0, l = 0)),
                                        axis.title.x = element_text(family = "Helvetica", size=22, margin = margin(t = 20, r = 0, b = 0, l = 0)),
                                        axis.title.y = element_text(family = "Helvetica", size=22, margin = margin(t = 0, r = 25, b = 0, l = 0)),
                                        axis.text = element_text(size=22)))
survdiff(Surv(AllKCur_rev$C_WaitTime, AllKCur_rev$Observed)~ AllKCur_rev$C_Environment, rho=0) #logrank test
summary(M1)
M1_se <- mean(summary(M1)$std.er)


## COX FORM: 
# UPDATE: the script for the model now reported in the paper (post revisions) is in extra_analysis_4reviews_eal
# ideal model: using coxme for random effects-- effect of environment & curiosity, random effects for ID & question number
All_SurvCur <- coxme(Surv(C_WaitTime, Observed) ~ C_Environment + CurZ + (1|ID) + (1|C_QuestionNum), data = AllKCur_rev)
All_SurvCur #(actual information to report)

## with environment & curiosity as predictors; USING THIS AS GRAPH FOR COXME MODEL
M2 <- coxph(Surv(C_WaitTime, Observed) ~ C_Environment + CurZ, data = AllKCur_rev)

M2_df <- with(AllKCur_rev,
              data.frame(C_Environment = c(0, 1), 
                         CurZ = rep(mean(CurZ, na.rm = TRUE), 2)
              )
)
M2_fit <- survfit(M2, newdata = M2_df)
ggsurvplot(M2_fit, data = M2_df, conf.int = TRUE, palette = c("#1f78b4", "#33a02c"), censor = FALSE, 
           legend = c(0.25, 0.15), legend.labs = c("Uniform Distribution", "Heavy-Tailed Distribution"), xlab="Wait Time (s)",
           ggtheme = theme_bw() + theme(legend.title = element_blank(), 
                                        legend.key.size = unit(1.5, 'lines'),
                                        legend.text=element_text(family = "Helvetica", size=18, margin = margin(t = 0, r = 25, b = 0, l = 0)),
                                        axis.title.x = element_text(family = "Helvetica", size=48, margin = margin(t = 20, r = 0, b = 0, l = 0)),
                                        axis.title.y = element_text(family = "Helvetica", size=48, margin = margin(t = 0, r = 25, b = 0, l = 0)),
                                        axis.text = element_text(size=30), panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank()))


########## AUC CURVES FOR SURVIVAL ANALYSIS
# Extracting AUC for each participant from the km_overall curves (Kaplan Meier curves)
print(km_ud_rev_overall, print.rmean=TRUE, rmean=20)$rmean
print(km_ht_rev_overall, print.rmean=TRUE, rmean=20)

print(M1, print.rmean=TRUE, rmean="common")

# Graphing AUC data
# load data 
setwd("/Users/emilylang/Dropbox/Research/ShohamyLab/KCur- Emily/Data")
AUC <- read.csv("KCur_JointAUC.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Do a paired t test
t.test(AUC ~ Environment, AUC, paired=TRUE)
se <- sd(AUC$Difference, na.rm=TRUE)/sqrt(64)

WS_Stat<- data.frame(Mean=c(16.14, 13.60), SE=c(se), Environment=c("UD", "HT"))
WS_Stat$Environment <- ifelse(WS_Stat$Environment == "UD", "Uniform Distribution", "Heavy-Tailed Distribution")
AUC$Environment <- ifelse(AUC$Environment == "UD", "Uniform Distribution", "Heavy-Tailed Distribution")

ggplot(data=WS_Stat, aes(x=Environment, y=Mean)) + 
  geom_point(data=AUC, aes(x=Environment, y=AUC, color=Environment, group=ID), size = 2) +
  geom_point(size=3) + geom_line(data=AUC, aes(x=Environment, y=AUC, group=ID), color= "dark grey") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, lwd=1.5) + 
  geom_line(aes(group=1), lwd=1) +
  scale_color_manual(values=c("#F8766D", "#00BFC4"), labels=c("Heavy-Tailed Distribution", "Uniform Distribution"), 
                     guide = guide_legend(reverse = TRUE)) + 
  labs(x="Environment", y="AUC (s)") + theme(legend.position= "none")


########## CURIOSITY ANALYSIS

## WILLINGNESS TO WAIT
## How does curiosity impact people's willingness to wait in either environment? (remove trials where wait
# is less than 1, but keep all other wait times >= 20)
invlog <- function(x) { 1/(1+exp(-(x))) }

KCur4 <- subset(AllKCur, AllKCur$C_WaitTime >= 1)

KCur4$C_Response <- recode(KCur4$C_Response, "S" = 0, "W" = 1)

#Exclude p's who waited for everything (temp, don't actually do this)
#KCur4 <- subset(KCur4, ID != 137 & ID != 129 & ID !=127 & ID != 125 & ID != 124 & ID !=112 & ID !=109 & ID !=102  & ID !=30  & ID !=20  & ID !=11)

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



####### MEMORY ANALYSIS
## Does curiosity predict memory? (excluding trials where wait was less than 1, but leaving all other
## wait times in); reported model!
invlog <- function(x) { 1/(1+exp(-(x))) }

M1 <- glmer(MC_Correct ~ CurZ * C_Environment + (CurZ + C_Environment |ID) + (C_Environment|C_QuestionNum),
            data=KCur4, family=binomial,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

summary(M1)

pred <- expand.grid(C_Environment = c(0, 1), CurZ=seq(min(KCur4$CurZ), max(KCur4$CurZ), by=.1))

fit <- as.data.frame(predict(M1, newdata=pred, re.form=NA))
fit <- cbind(pred, fit)
colnames(fit)[3] <- "fit"
fit$Prob <- invlog(fit$fit)

Condse <- summary(M1)$coefficients[2,2]
fit$se <- Condse
fit$lower <- fit$fit - fit$se*2
fit$upper <- fit$fit + fit$se*2
fit$lowerprob <- invlog(fit$lower)
fit$upperprob <- invlog(fit$upper) 

ggplot(fit, aes(x=CurZ, y=Prob, color=as.factor(C_Environment), group=as.factor(C_Environment))) + 
  geom_ribbon(aes(ymin=lowerprob, ymax=upperprob, fill = as.factor(C_Environment)), alpha=.3, color=NA) +
  geom_line(lwd=1.5) + scale_color_manual(values=UD_HT_colors, labels=c("Uniform Distribution", "Heavy-Tailed Distribution")) +
  scale_y_continuous(lim=c(0.4,1)) + labs(x="Question Curiosity (Z)", y="P(Remember)") +
  scale_fill_manual(values = UD_HT_colors, labels = c("Uniform Distribution", "Heavy-Tailed Distribution")) 

## Do curiosity and wait time interact to predict memory?
##(excluding trials where wait was less than 1, but leaving all other wait times in)
## Use wait time rather than environment as a predictor

M3 <- glmer(MC_Correct ~ CurZ * C_WaitTime + (CurZ * C_WaitTime |ID),
            data=KCur4,
            family=binomial,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

summary(M3)

#model with median split to plot continuous interaction
KCur4$WaitTime <- 0
median <- median(KCur4$C_WaitTime)
KCur4$WaitTime[which(KCur4$C_WaitTime > median)] <- 1

M4 <- glmer(MC_Correct ~ CurZ * WaitTime + (CurZ * WaitTime |ID),
            data=KCur4,
            family=binomial,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

summary(M4)

pred4 <- expand.grid(WaitTime = c(0, 1), CurZ=seq(min(KCur4$CurZ), max(KCur4$CurZ), by=.1))

fit4 <- as.data.frame(predict(M4, newdata=pred4, re.form=NA))
fit4 <- cbind(pred4, fit4)
colnames(fit4)[3] <- "fit4"
fit4$Prob <- invlog(fit4$fit4)

Condse4 <- summary(M4)$coefficients[2,2]
fit4$se <- Condse4
fit4$lower <- fit4$fit4 - fit4$se*2
fit4$upper <- fit4$fit4 + fit4$se*2
fit4$lowerprob <- invlog(fit4$lower)
fit4$upperprob <- invlog(fit4$upper)

ggplot(fit4, aes(x=CurZ, y=Prob, color=as.factor(WaitTime), group=as.factor(WaitTime))) + 
  geom_ribbon(aes(ymin=lowerprob, ymax=upperprob, fill = as.factor(WaitTime)), alpha=.3, color=NA) +
  geom_line(lwd=1.5) + scale_color_manual(values=UD_HT_colors, labels=c("Long Wait", "Short Wait")) +
  scale_y_continuous(lim=c(0.4,1)) + labs(x="Question Curiosity (Z)", y="P(Remember)") +
  scale_fill_manual(values = UD_HT_colors, labels = c("Long Wait", "Short Wait")) 

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


