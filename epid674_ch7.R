################################################
# EPID674 Epidemiologic Data Analysis using R 
# R Script for Chapter 7. Other Useful Packages
################################################

##### Load these packages all the time
library(foreign)
library(Hmisc)
library(epicalc)

setwd("C:/Desktop/epid798")	# set working directory
#setwd("/Users/kbakulsk/Google Drive/Teaching/EPID798/Lab") #my example 2: Google Drive on my laptop
#setwd("M:/EPID798/Lab") #my example 3: Umich M-drive
#setwd("/afs/umich.edu/user/b/a/bakulski/EPID798") #my example 4: Umich AFS storage

# Load the save R data and attach it
load("nhanes3_v2.rda")
attach(nhanes3)

##########################################################
## 7.1. nlme: linear and nonlinear mixed effects models ##
##########################################################

## Try linear mixed effects models using nlme
library(nlme)

##Distribution of primary sampling unit (psu) and cluster (strata)
tab1(nhanes3$psu)		# 2 PSU
tab1(nhanes3$strata)	# 49 strata

nhanes3$locode<-ifelse(nhanes3$psu==1, nhanes3$strata, nhanes3$strata+49)
tab1(nhanes3$locode)
table(nhanes3$psu,nhanes3$locode)

###Random intercepts only model for SBP
sbp.lme<-lme(fixed=log(sbp)~log(bpb)+age+I(age^2)+bmi+factor(race)+factor(sex)
+factor(smk)+factor(educ)+factor(alc), random=~1|locode, na.action=na.omit, data=nhanes3)
summary(sbp.lme)

##Return variance-covariance matrix for random variables
getVarCov(sbp.lme)

##Random Intercepts and slopes model
sbp.lme2<-lme(fixed=log(sbp)~log(bpb)+age+I(age^2)+bmi+factor(race)+factor(sex)
+factor(smk)+factor(educ)+factor(alc), random=~1+log(bpb)|locode, na.action=na.omit, data=nhanes3)
summary(sbp.lme2)
getVarCov(sbp.lme2)

## compare these two models
anova(sbp.lme, sbp.lme2)

#########################################################
#### 7.2. survey: analysis of complex survey samples ####
#########################################################

#### load the survey package
#install.packages("survey")
library(survey)
library(splines)

#### first specify the design effect
bpdsn<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=nhanes3, nest=T)
#nest=T: relabel cluster ids to enforce nesting within strata

#### Descriptive analyses
svymean(~age, design=bpdsn)
summ(age)
svymean(~age+sbp+bpb, design=bpdsn)
svymean(~sbp+bpb+bmi, design=bpdsn)
svymean(~sbp+bpb+bmi, design=bpdsn, na.rm=T)

svytable(~sex, design=bpdsn)
svytable(~sex, Ntotal=100, design=bpdsn)
svytable(~sex+race, design=bpdsn)
summary(svytable(~sex+race, Ntotal=100, design=bpdsn))
svychisq(~sex+race, design=bpdsn)
chisq.test(sex, race)
svytable(~htn+smk, design=bpdsn)
svychisq(~htn+smk, design=bpdsn, statistic="Wald")
## Survey Chisq used in SUDAAN

## Univariate distributions and bivariate associations using graphical analyses
svyhist(~sbp, design=bpdsn)
svyhist(~log(sbp), design=bpdsn)
svyhist(~bpb, bpdsn)
svyhist(~log(bpb), bpdsn)

svyboxplot(log(sbp)~1, bpdsn)
svyboxplot(log(sbp)~as.factor(smk), bpdsn)

bp.bysexrace<-svyby(~sbp+dbp, ~sex+race, bpdsn, svymean)
barplot(bp.bysexrace)
barplot(bp.bysexrace, legend=TRUE)

## change the label of x-axis and ylim
xlabel<-c("Male, White","Female, White","Male, Black","Female, Black")
barplot(bp.bysexrace, legend=TRUE, ylim=c(0,160),col=c("purple","violet"),
names=xlabel)

#simple scatterplot with circles whose area is proportional to the sampling weight
svyplot(log(sbp)~age, bpdsn, style="bubble")
#style="transparent" plots points with opacity proportional to sampling weight 
svyplot(log(sbp)~log(bpb), bpdsn, style="transparent", pch=19, xlab="Blood lead", ylab="log(SBP)")

## Scatterplot smoothing
smth.sbp.bpb1<-svysmooth(log(sbp)~bpb, bpdsn, bandwidth=10)
#fit local polynomial kernel smoothing with a bandwidth=10 bpb unit 
plot(smth.sbp.bpb1)
smth.sbp.bpb2<-svysmooth(log(sbp)~bpb, bpdsn, bandwidth=20)
plot(smth.sbp.bpb2)


#### Linear and mixed effect models
## First fit a linear regression
sbp.lm<-lm(log(sbp)~log(bpb)+age+I(age^2)+bmi+factor(race)
+factor(sex)+factor(smk)+factor(educ)+factor(alc),data=nhanes3)
summary(sbp.lm)

sbp.lme<-lme(fixed=log(sbp)~log(bpb)+age+I(age^2)+bmi+factor(race)
+factor(sex)+factor(smk)+factor(educ)+factor(alc), random=~1|strata, 
na.action=na.omit, data=nhanes3)
summary(sbp.lme)

## Let's use the survey package
sbp.svy<-svyglm(log(sbp)~log(bpb)+age+I(age^2)+bmi+factor(race)
+factor(sex)+factor(smk)+factor(educ)+factor(alc),bpdsn)
summary(sbp.svy)
plot(sbp.svy)

par(mfrow=c(2,2))
plot(sbp.lm, which=1)
plot(sbp.lm, which=2)
plot(sbp.svy, which=1)
plot(sbp.svy, which=2)
par(mfrow=c(1,1))

## compare the beta's
summary(sbp.lm)$coef[2,]
summary(sbp.lme)$tTable[2,]
summary(sbp.svy)$coef[2,]

#### Logistic regression model
htn.svy<-svyglm(htn~log(bpb)+ns(age,df=5)+ns(bmi,df=5)+factor(race)
+factor(sex)+factor(educ)+hematoc+chol+packyrs+diag_dm,
family=quasibinomial(),bpdsn)
summary(htn.svy)

########################################
#### 7.3. rmeta: Meta-analysis in R ####
########################################

#install.packages("rmeta")
library(rmeta)

##Suppose you have the estimates from 16 cities of the effect of PM2.5 
##on Myocardial Infarction (MI) hospital admissions. 
##We want to obtain the combined effect across all the cities.
city<-1:16
coef<-c(0.00155,0.00445,-0.00597,0.00237,0.00031,-0.00035,0.00206,0.00127,
-0.00395,-0.00434,0.00241,0.00730,-0.00175,0.00117,0.00007,0.00350)

se<-c(0.013589379,0.003224905,0.002400626,0.001797566,0.001019834,0.002658267,
0.001529252,0.002388748,0.002583907,0.005885021,0.00122501,0.00597465,
0.002546777,0.002406304,0.000913846,0.00517299)

model<-meta.summaries(coef, se, names=city, method="random")
#method="random" estimates and adds a heterogeneity variance
summary(model)

##Plot the city-specific and summary results
par(mfrow=c(1,1))
plot(model)
metaplot(coef, se)

##Get combined estimate, its standard error and RR and 95% CI
model$summary 
model$summ
model$se.summary 
model$se

exp(model$summ*10)
exp((model$summ-1.96*model$se)*10)
exp((model$summ+1.96*model$se)*10)

### End of the class ###
### Thank you so much! ###
