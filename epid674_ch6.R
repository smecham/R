###############################################
# EPID674 Epidemiologic Data Analysis using R 
# R Script for Chapter 6. Regression Models
###############################################

##### Load these packages all the time
library(foreign)
library(Hmisc)
library(epicalc)
library(stargazer)

setwd("C:/Desktop/epid798")	# set working directory
#setwd("/Users/kbakulsk/Google Drive/Teaching/EPID798/Lab") #my example 2: Google Drive on my laptop
#setwd("M:/EPID798/Lab") #my example 3: Umich M-drive
#setwd("/afs/umich.edu/user/b/a/bakulski/EPID798") #my example 4: Umich AFS storage

### You can start R afresh if you didn't save the workspace image
# Load the save R data
load("nhanes3.rda")

# attach our data
attach(nhanes3)

### 6.1. Linear Models: Association between systolic blood pressure (SBP) and blood lead (bpb)

##Does the distribution of log(sbp) look closer to the normal distribution?
par(mfrow=c(2,1))
hist(sbp,nclass=20, col="darkorchid")
hist(log(sbp),nclass=20, col="seagreen3")

## Let?s start with non-log transformed SBP first. Look at bivariate association between sbp and continuous covariates
par(mfrow=c(1,1))
plot(age,sbp, pch=20, cex=0.7, col="dodgerblue",cex.lab=1.2, las=1, ylab="SBP", xlab="Age (years)")
lines(smooth.spline(age,sbp, df=10), col = "dodgerblue", lwd=3)
plot(bmi,sbp, pch=20, cex=0.7, col="dodgerblue",cex.lab=1.2, las=1, ylab="SBP", xlab="BMI (kg/m2)")
lines(smooth.spline(na.omit(bmi),sbp[na.omit(bmi)], df=3, ), col = "grey60", lwd=3)
plot(bpb,sbp, pch=20, cex=0.7, col="dodgerblue",cex.lab=1.2, las=1, ylab="SBP", xlab="Blood lead level (ug/dL)")
lines(smooth.spline(bpb,sbp, df=10), col = "grey60", lwd=3)

#Clean up dataset. Make sure NaN are NA for key covariates
table(nhanes3$educ, useNA = "always")
nhanes3$educ[is.nan(nhanes3$educ)]<-NA
table(nhanes3$educ, useNA = "always")
table(nhanes3$alc, useNA = "always")
nhanes3$alc[is.nan(nhanes3$alc)]<-NA
table(nhanes3$alc, useNA = "always")
attach(nhanes3)


##Start creating a simple regression model for systolic blood pressure, 
##including only blood lead (bpb) (crude model).
sbp.model<-lm(sbp~bpb, na.action=na.omit, data=nhanes3)
summary(sbp.model)
summary.aov(sbp.model)
anova(sbp.model)

##Add age in the model
sbp.model1<-lm(sbp~bpb+age, na.action=na.omit, data=nhanes3)
summary(sbp.model1)

##Add race, which is a categorical variable, therefore we should use factor()
##or as.factor() which creates indicator variables for each category of race
table(race)

sbp.model2<-lm(sbp~bpb+age+factor(race), na.action=na.omit, data=nhanes3)
summary(sbp.model2)

####Change the reference level for a factor variable
#Method 1: Make a new (or alter the old) variable
table(race)
race2<-relevel(factor(race), ref=2)
table(race2)
sbp.model2<-lm(sbp~bpb+age+factor(race2), na.action=na.omit, data=nhanes3)
summary(sbp.model2)
#Method 2: Change the contrasts in the lm() statement
sbp.model2<-lm(sbp~bpb+age+C(factor(race),contr.treatment(2, base=2)), na.action=na.omit, data=nhanes3)
summary(sbp.model2)

# R provides Type I sequential SS, not the default Type III marginal SS reported by SAS and SPSS. 
# We will need to use the drop1() function to produce the familiar Type III results. 
# It will compare each term with the full model. 
anova(sbp.model2)
drop1(sbp.model2, test="F")

## alternative
library(car)
Anova(sbp.model2, type="III")

##We can also use the update function to add variables to a pre-existing model
sbp.model2<-update(sbp.model1,.~.+ factor(race))
summary(sbp.model2)

##Add other covariates to the model: which variables are biologically important?
##Add sex, BMI, educ and smk.
sbp.model3<-lm(sbp~bpb+age+factor(race)+factor(sex)+bmi+factor(educ)+factor(smk), na.action=na.omit, data=nhanes3)
summary(sbp.model3)

#See output from models 1,2,and 3 side by side
stargazer(sbp.model1, sbp.model2, sbp.model3, type="text", dep.var.labels = "Systolic Blood Pressure (mmHg)")

### Check if alcohol consumption is a confounder
sbp.model4<-update(sbp.model3,.~.+factor(alc))
summary(sbp.model4)

summary(sbp.model4)$coef
summary(sbp.model4)$coef[2,1]
(summary(sbp.model4)$coef[2,1]-summary(sbp.model3)$coef[2,1])/summary(sbp.model3)$coef[2,1]


## Compute percent increase in SBP and 95% CI per one unit increase and IQR increase in bpb
# for one unit increase
change<-summary(sbp.model4)$coef[2,1]
l95ci<-summary(sbp.model4)$coef[2,1]-1.96*summary(sbp.model4)$coef[2,2]
u95ci<-summary(sbp.model4)$coef[2,1]+1.96*summary(sbp.model4)$coef[2,2]
change
l95ci
u95ci
regress.display(sbp.model4)

# for an IQR increase
IQR(bpb)
change.iqr<-IQR(bpb)*summary(sbp.model4)$coef[2,1]
l95ci.iqr<-IQR(bpb)*(summary(sbp.model4)$coef[2,1]-1.96*summary(sbp.model4)$coef[2,2])
u95ci.iqr<-IQR(bpb)*(summary(sbp.model4)$coef[2,1]+1.96*summary(sbp.model4)$coef[2,2])
change.iqr
l95ci.iqr
u95ci.iqr
regress.display(sbp.model4)$table[2,]*IQR(bpb)


##Does the variable packyrs give a better estimate than smk?
sbp.model5<-update(sbp.model4,.~.-factor(smk)+packyrs)
summary(sbp.model5)

AIC(sbp.model4, sbp.model5)

##Run two models one with age and one adding age as a quadratic function. 
##Does the quadratic term improve the fit of the model?
sbp.model6<-update(sbp.model4,.~.+I(age^2))
summary(sbp.model6)

##use the anova function to compare the 2 models and 
##see if the quadratic term improves the model.
anova(sbp.model4, sbp.model6, test="F")

## so sbp.model6 is our final model
## What is the relationship between bpb and sbp?


###Regression Diagnostics
##In the case of linear model, the plot of the model gives diagnostic plots
par(mfrow=c(1,1))
plot(sbp.model6)

par(mfrow=c(2,2))
plot(sbp.model6, id.n = 5, cex=0.1)
#plot.lm(sbp.model6)

par(mfrow=c(1,1))
plot(sbp.model6, which=1)
plot(sbp.model6, which=2)
plot(sbp.model6, which=3)
plot(sbp.model6, which=5)

plot(sbp.model6, which=4)

## How about log(sbp)?
## Check how diagnostic plots using log-transformed sbp look like

sbp.model6.log<-update(sbp.model6,.-sbp+log(sbp)~.)
summary(sbp.model6.log)
plot(sbp.model6.log)
hist(residuals(sbp.model6.log),nclass=20)

par(mfrow=c(2,2))
plot(sbp.model6, which=1)
plot(sbp.model6.log, which=1)
hist(residuals(sbp.model6), main="Histogram of res(SBP)")
hist(residuals(sbp.model6.log), main="Histogram of res(log(SBP))")
par(mfrow=c(1,1))

###Partial Residual Plots
##Plot sbp vs. bpb given that other variables are in the model (adjusted)

##This can be done by 'termplot'
termplot(sbp.model6, partial.resid=TRUE, col.res="gray30", cex=0.3,  smooth=panel.smooth)


## Suppose you decide to go with log transformed SBP. 

## Compute percent increase in SBP and 95% CI per one unit increase and IQR increase in bpb
## what to report if y is log-transformed?
## percent increase (difference) vs. absolute increase (difference)?
## Compute percent increase in SBP and 95% CI per one unit increase and IQR increase in bpb

summary(sbp.model6.log)$coef
summary(sbp.model6.log)$coef[2,1]
summary(sbp.model6.log)$coef[2,2]
exp(summary(sbp.model6.log)$coef[2,1])
#exp(sbp.model6.log$coef[2])
p.change<-100*(exp(summary(sbp.model6.log)$coef[2,1])-1)
l95ci<-100*(exp(summary(sbp.model6.log)$coef[2,1]-1.96*summary(sbp.model6.log)$coef[2,2])-1)
u95ci<-100*(exp(summary(sbp.model6.log)$coef[2,1]+ 1.96*summary(sbp.model6.log)$coef[2,2])-1)
p.change
l95ci
u95ci

IQR(bpb)
p.change.iqr<-100*(exp(IQR(bpb)*summary(sbp.model6.log)$coef[2,1])-1)
l95ci.iqr<-100*(exp(IQR(bpb)*(summary(sbp.model6.log)$coef[2,1]-1.96*summary(sbp.model6.log)$coef[2,2]))-1)
u95ci.iqr<-100*(exp(IQR(bpb)*(summary(sbp.model6.log)$coef[2,1]+ 1.96*summary(sbp.model6.log)$coef[2,2]))-1)
p.change.iqr
l95ci.iqr
u95ci.iqr


## linearity assumption
## smoothing plot

library(mgcv)

sbp.model6.gam<-gam(sbp~s(bpb)+age+factor(race)+factor(sex)+bmi+factor(educ)+factor(smk)+factor(alc)+I(age^2), na.action=na.omit, data=nhanes3)
summary(sbp.model6.gam)
plot(sbp.model6.gam)

plot(sbp.model6.gam, xlab="blood lead (ug/dL)", ylab="Change in log(SBP)")

sbp.model7<-lm(sbp~log(bpb)+age+factor(race)+factor(sex)+bmi+factor(educ)+factor(smk)+factor(alc)+I(age^2), na.action=na.omit, data=nhanes3)
summary(sbp.model7)
summary(sbp.model6)

## compare model6 and model7
AIC(sbp.model6, sbp.model7)


###Effect modification by gender
table(sex)

### Stratified by gender
##Male (sex==1)
sbp.model6.male<-lm(sbp~bpb+age+factor(race)+bmi+factor(educ)
               +factor(smk)+factor(alc)+I(age^2),data=nhanes3, subset=(sex==1))
summary(sbp.model6.male)

##Female (sex==2)
sbp.model6.female<-lm(sbp~bpb+age+factor(race)+bmi+factor(educ)
                    +factor(smk)+factor(alc)+I(age^2),data=nhanes3, subset=(sex==2))
summary(sbp.model6.female)

### Use interaction model
sbp.model6.int<-lm(sbp~bpb+factor(sex)+bpb*factor(sex)+age+factor(race)+bmi+factor(educ)
                   +factor(smk)+factor(alc)+I(age^2),data=nhanes3)
# below is the same
sbp.model6.int<-lm(sbp~bpb*factor(sex)+age+factor(race)+bmi+factor(educ)
                   +factor(smk)+factor(alc)+I(age^2),data=nhanes3)
summary(sbp.model6.int)

########################################
############ Exercise 6.2 ##############
########################################

####################################################################################################

### 6.3. Generalized Linear Models: Association between hypertension (htn) and blood lead (bpb)

##Logistic regression for hypertension
##Look at hypertension (htn)
tab1(htn, graph=F)

htn.model<-glm(htn~bpb+age+factor(sex)+factor(race)+bmi+factor(educ)+factor(smk)+factor(alc), family=binomial, 
               na.action=na.omit, data=nhanes3)
summary(htn.model)

## you can construct your model using variable selection methods (forward, backward, stepwise)
# first, define a matrix of predictors
X<-nhanes3[,c(3:8, 13:25, 29, 30)]
Y<-htn
dat=data.frame(cbind(Y,X))
dat<-na.omit(dat)	#note: dataset should be complete
summ(dat)

fit.start=glm(Y~1,data=dat, family=binomial)
summary(fit.start)
fit.full=glm(Y~.,data=dat, family=binomial)
summary(fit.full)

# forward
full<-formula(glm(Y~.,data=dat, family=binomial))
full
fit.forward=step(fit.start,direction='forward',scope=full)
summary(fit.forward)

# if you want to keep bpb
fit.forward1=step(fit.start,direction='forward',scope=list(lower=~bpb, upper=full))
fit.forward1=step(glm(Y~bpb, family=binomial, data=dat),direction='forward', scope=list(upper=full))
summary(fit.forward1)

# backward (always keep bpb in the model)
fit.backward=step(fit.full,direction='backward', scope=list(lower=~bpb))
summary(fit.backward)

# stepwise (always keep bpb in the model)
fit.step=step(fit.full,direction='both', scope=list(lower=~bpb))
summary(fit.step)		# this is exactly the same as backward

fit.step1=step(glm(Y~bpb, family=binomial, data=dat),direction='both', scope=list(lower=~bpb, upper=full))
summary(fit.step1)


# compute ORs
logistic.display(htn.model)

# regression diagnostic
plot(htn.model)
plot(htn.model, which=4)

par(mfrow=c(2,2))
plot(htn.model)

par(mfrow=c(1,1))
termplot(htn.model)
termplot(htn.model, se=T)
termplot(htn.model, se=T, partial.resid=T)

####################################################################################################
##Poisson regression for respiratory death: Montana from epicalc
data(Montana)
summ(Montana)
head(Montana, 10)
hist(Montana$respdeath)

par(mfrow=c(2,2))
tab1(Montana$agegr)
tab1(Montana$period)
tab1(Montana$start)
tab1(Montana$arsenic)

# label categorical variables
Montana$agegr<-factor(Montana$agegr, labels=c("40-49","50-59","60-69","70-79"))
Montana$period<-factor(Montana$period, labels=c("1938-1949", "1950-1959", "1960-1969", "1970-1977"))
Montana$start<-factor(Montana$start, labels=c("pre-1925", "1925 & after"))
Montana$arsenic<-factor(Montana$arsenic, labels=c("<1 year", "1-4 years","5-14 years", "15+ years"))

tab1(Montana$agegr, missing=F)
tab1(Montana$period, missing=F)
tab1(Montana$start, missing=F)
tab1(Montana$arsenic, missing=F)
par(mfrow=c(1,1))

# Compute incidence rate by age and period
table.pyears<-tapply(Montana$personyrs, list(Montana$period, Montana$agegr), sum)
table.deaths<-tapply(Montana$respdeath, list(Montana$period, Montana$agegr), sum)
table.inc10000<-table.deaths/table.pyears*10000
table.inc10000

# create a time-series plot of the incidence
plot.ts(table.inc10000, plot.type="single", xlab="", ylab="#/10,000 person-years", xaxt="n", col=c("black","blue","red","green"), lty=c(2,1,1,2), las=1)
points(rep(1:4,4), table.inc10000, pch=22, cex=table.pyears/sum(table.pyears)*20)
title(main = "Incidence by age and period")
axis(side = 1, at = 1:4, labels = levels(Montana$period))
legend("topleft", legend=levels(Montana$agegr)[4:1], col=c("green","red", "blue", "black"),bg="white",lty=c(2,1,1,2))

# check arsenic
tab1(Montana$arsenic)
tapply(Montana$respdeath, Montana$arsenic, mean)
tapply(Montana$personyrs, Montana$arsenic, mean)

#Poisson model
resp.mode11<-glm(respdeath~period, offset=log(personyrs), family=poisson, data=Montana)
summary(resp.mode11)

resp.mode12<-glm(respdeath~agegr, offset=log(personyrs), family=poisson, data=Montana)
summary(resp.mode12)

resp.mode13<-glm(respdeath~period+agegr, offset=log(personyrs), family=poisson, data=Montana)
summary(resp.mode13)

AIC(resp.mode11, resp.mode12, resp.mode13)
## model2 is better

resp.mode14<-glm(respdeath~agegr+arsenic, offset=log(personyrs), family=poisson, data=Montana)
summary(resp.mode14)

# is there a linear trend across arsenic exposure?
resp.mode14.lin<-glm(respdeath~agegr+as.numeric(arsenic), offset=log(personyrs), family=poisson, data=Montana)
summary(resp.mode14.lin)

## compute IRR
idr.display(resp.mode14)


########################################
############ Exercise 6.4 ##############
########################################


####################################################################################################

### 6.5. Matched Case-Control Study: VC1to1 from epicalc

data(VC1to1)
summ(VC1to1)
head(VC1to1)

# Reshape the data to facilitate data exploration
# function 'reshape' converts wide to long or vice versa

wide <- reshape(VC1to1, timevar="case", v.names=c("smoking","rubber", "alcohol"), idvar="matset", direction="wide")
head(wide,3)
table(wide$smoking.1, wide$smoking.0, dnn=c("smoking in case", "smoking in control"))
#dnn: dimnames names 

# matchTab() is for the conditional OR (McNemar's OR)

matchTab(VC1to1$case, VC1to1$smoking, strata=VC1to1$matset)
matchTab(VC1to1$case, VC1to1$rubber, strata=VC1to1$matset)
matchTab(VC1to1$case, VC1to1$alcohol, strata=VC1to1$matset)

## look at the full dataset VC1to6

data(VC1to6)
summ(VC1to6)
VC1to6[,]

# what is the effect of smoking?
matchTab(VC1to6$case, VC1to6$smoking, strata=VC1to6$matset)
matchTab(VC1to6$case, VC1to6$alcohol, strata=VC1to6$matset)

# conditional logistic reg using clogit from survival package

library(survival)
clogit1<-clogit(case~smoking+alcohol+strata(matset), data= VC1to1)
summary(clogit1)

clogit2<-clogit(case~smoking+alcohol+strata(matset), data= VC1to6)
summary(clogit2)

# compute ORs
clogistic.display(clogit1)
clogistic.display(clogit2)


####################################################################################################

### 6.6. Survival Analysis: Association between total mortality (d_total) and blood lead (bpb)

### load the survival package
library(survival)

tab1(d_total)
summ(pmon_mec)

### Define Surv()
surv.total<-Surv(pmon_mec, d_total)
surv.total

### K-M Life table and curve
fit.total<-survfit(Surv(pmon_mec, d_total)~1)
summary(fit.total)

plot(fit.total)
## suppress 95% CI lines and the time marks for censored subjects.
plot(fit.total, ylim=c(0.7,1.0), conf.int=F, mark.time=F)

### Survival by different levels of covariates
fit.total.sex<-survfit(Surv(pmon_mec, d_total)~sex)
fit.total.sex
summary(fit.total.sex)

plot(fit.total.sex, col=c("blue","red"), lty=c(1,2))
plot(fit.total.sex, ylim=c(0.6,1.0), col=c("blue","red"), lty=c(1,2), mark.time=F)
title(main="Kaplan-Meier curve", xlab="Time (months)", ylab="Survival probability")
legend("topright", legend=c("Men","Women"), lty=c(1,2), col=c("blue","red"))

### Test for differences in survival curves
survdiff(Surv(pmon_mec, d_total)~sex)

### Cox regression

cox.bpb<-coxph(Surv(pmon_mec, d_total)~bpb)
summary(cox.bpb)

bpb3<-cut2(bpb, g=3)
tab1(bpb3)

# K-M Life table and curve
fit.total.bpb3<-survfit(Surv(pmon_mec, d_total)~bpb3)
summary(fit.total.bpb3)
plot(fit.total.bpb3, col=c(1:3), lty=c(1:3), mark.time=F)
plot(fit.total.bpb3, col=c(1:3), lty=c(1:3), mark.time=F, ylim=c(0.6,1.0))
title(main="Survival curve in relation to blood lead levels", xlab="Time (months)", ylab="Survival probability")
legend(30,0.7, legend=c("Q1","Q2","Q3"), lty=c(1:3), col=c(1:3))

#crude
cox.bpb3<-coxph(Surv(pmon_mec, d_total)~bpb3)
summary(cox.bpb3)

#adjusted
cox.bpb3.adj<-coxph(Surv(pmon_mec, d_total)~bpb3+age+factor(sex)+factor(race)+factor(educ)+factor(smk)+factor(alc))
summary(cox.bpb3.adj)

### Test for the proportional hazards assumption
test.prop<-cox.zph(cox.bpb3.adj)
test.prop

## Display a graph of the scaled Schoenfeld residuals, along with a smooth curve
plot(test.prop)  # for all variables
plot(test.prop, var=1)
plot(test.prop, var=2)
plot(test.prop, var=4)
abline(h=0, lty=3, col=2)

# stratify by sex
cox.bpb3.adj1<-coxph(Surv(pmon_mec, d_total)~bpb3+age+strata(sex)+factor(race)+factor(educ)+factor(smk)+factor(alc))
summary(cox.bpb3.adj1)

test.prop1<-cox.zph(cox.bpb3.adj1)
test.prop1

####################################################################################################
### 6.7. Write your own functions

## Let's make a simple macro that calculates the mean and standard deviation at the same time.

mystats<-function(x)
{
  print(mean(x, na.rm=T))
  print(sd(x, na.rm=T))
}

mystats(age)
mystats(bpb)

# get the results in vector form
mystats<-function(x)
{
  mymean<-mean(x, na.rm=T)
  mysd<-sd(x, na.rm=T)
  c(mean=mymean, sd=mysd)
}

mystats(age)
mystats(bpb)	


# Assume that you are examining age-adjusted associations of SBP with from the 13th variables (bmi) to 25th variables (packyrs) (n=13)
# you want to run 13 linear regression models and save beta's and p-values

summ(nhanes3)
test.var<-nhanes3[,c(13:25)]
head(test.var)

# example
mod<-lm(sbp~bmi+age, data=nhanes3, na.action=na.omit)
summary(mod)

# how to extract beta for bmi?
summary(mod)$coef[2,1]

# how to extract p-value for bmi?
summary(mod)$coef[2,4]

Test<-function(data, y, cov){
  nvar<-ncol(data)
  newdata<-data.frame(cbind(data, y, cov))
  
  tmatrix<-data.frame(matrix(NA,2,nvar)) # 2 rows  
  colnames(tmatrix)<-colnames(data)  # create a row for each var  
  rownames(tmatrix)<-c("beta","p")   
  
  for(i in 1:nvar){ 
    ind<-data[,i] 
    model<-lm(y~ind+cov, data=newdata, na.action=na.omit) 
    tmatrix[1,i]<- summary(model)$coef[2,1] 
    tmatrix[2,i]<- summary(model)$coef[2,4] 
  } 
  return(tmatrix) 
  write.csv(tmatrix, file='sbp.results.csv') 
}  

Test(test.var, sbp, age)


