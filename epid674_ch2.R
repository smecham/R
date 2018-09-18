###############################################
# EPID674 Epidemiologic Data Analysis using R 
# R Script for Chapter 2, Exploring data with R
###############################################

### Comments can be put behind the # symbol 
### Anything from the # to the end of the line will be ignored by R  

sessionInfo() #check what packages are already loaded by default

# Install packages. Do this only once.
install.packages("Hmisc")#if a package has trouble loading, you will need to install it. Do this only once
install.packages("epitools")
install.packages("Epi")
install.packages("sas7bdat") 
install.packages("gmodels")
install.packages("abind")
install.packages("stargazer")
install.packages("https://cran.r-project.org/src/contrib/Archive/epicalc/epicalc_2.15.1.0.tar.gz", repos=NULL)


# Load packages. Load relevant packages every time you start a new R session
library(foreign)
library(Hmisc)
library(epitools)
library(Epi)
library(sas7bdat)
library(abind)
library(gmodels)
library(epicalc)
library(stargazer)

# Set your working directory (default location for loading/saving files)
setwd("C:/Desktop/epid798")	# Change location to fit your computer: Folder on your computer
#setwd("/Users/kbakulsk/Google Drive/Teaching/EPID798/Lab") #my example 2: Google Drive on my laptop
#setwd("M:/EPID798/Lab") #my example 3: Umich M-drive
#setwd("/afs/umich.edu/user/b/a/bakulski/EPID798") #my example 4: Umich AFS storage
getwd() #check your working directory

# Load data from SAS
nhanes3<-read.sas7bdat("nhanes3.sas7bdat")
colnames(nhanes3)

# Save data as an R object
save(nhanes3, file="nhanes3.rda")
load("nhanes3.rda") # Try reading in the data that you just saved

# Exporting data as a txt/csv file
write.table(nhanes3, file="nhanes3.txt")
write.csv(nhanes3, file="nhanes3.csv")

# Try reading in the .csv that you just made
nhanes3c<-read.csv("nhanes3.csv")
names(nhanes3c)

## list all objects and remove the ones we do not want
ls()
rm(x,y,z) 
#rm(list=ls())  #removes all objects
ls()

# Explore the data set
class(nhanes3) #What type of object is it?
dim(nhanes3) #What are the dimensions?
names(nhanes3) #What are the column names?
acorn(nhanes3) #What do you see here?
head(nhanes3) #What do the first 6 rows look like?
nhanes3[1:10,] #What do the first 10 rows look like?

# Explore the variables
length(nhanes3$age)

attach(nhanes3) #attaches the dataset so you can use the column names without specifying the data.frame

### Create a factor variable (convert from a numeric variable)
class(sex)
sex1 <- as.factor(sex)	## creating sex1 from sex
class(sex1)
is.factor(sex)
is.factor(sex1)
levels(sex1) <- c("male", "female")
sex[1:10]
sex1[1:10]

# OR
sex1 <- factor(sex, levels=c(1,2), labels=c("male", "female"))
sex1 <- factor(sex, labels=c("male", "female"))
sex1[1:10]
table(sex1)
table(sex)

### opposite way is as.numeric()
sex2<-as.numeric(sex1)
sex2[1:10]
table(sex2)

### Data description functions
mean(x=age, trim=0, na.rm=FALSE)
mean(na.rm=FALSE, x=age, trim=0)
mean(age, 0, F)
mean(age)

median(age)
quantile(age)
quantile(age, c(0.1,0.9))
sd(age)
IQR(age)

summary(age)
summary(nhanes3)

########################################
############ Do Exercise 2A ##############
########################################

## Calculate descriptives on a subset of the dataset
summary(nhanes3[,1:5])
summary(nhanes3[,16:20])
summary(nhanes3[,c(1,6,11,16,21)])
summary(nhanes3[,c("age","sex","race","bpb","sbp")])
summary(nhanes3[,-c(1:30)])

## Identify the positions of each TRUE values where the condition is true using which()
which(age>80)
which(bmi<18)

# Using the which() function,  you can also run any R commands for subsets.
table(race)
white<-which(race==1)
length(white)
black<-which(race==2)
length(black)

malewhite<-which(nhanes3$sex==1&nhanes3$race==1)
length(malewhite)

# subset(): You can create a subset from a data frame using the subset() function. 
# For example, create a dataset for race==1 (White).
nhanes3.w<-subset(nhanes3, race==1)
summary(nhanes3.w)

# To select only the 'age', 'bmi', and 'sex' variables among White people
nhanes3.w<-subset(nhanes3, race==1, select=c(age,bmi,sex))
summary(nhanes3.w)

## sample() to randomly sample
samp<-sample(seqn, 100)
samp
samp1<-sample(seqn, 100, replace=T)
samp1

## create a categorical variable
AGE5a<-cut(age, 5)	# five equally space intervals (14 yrs)
summary(AGE5a)
AGE5b<-cut(age, quantile(age, c(0,.2,.4,.6,.8,1)), include.lowest=T) # quintiles
summary(AGE5b)
AGE5c<-cut(age, breaks=c(19,40,50,60,70,90))
summary(AGE5c)

AGE5d<-cut2(age, g=5)
summary(AGE5d)
AGE5e<-cut2(age, c(40,50,60,70)) 
summary(AGE5e)

#Another way
agecat <- age
agecat[age<30] <- "20's"
agecat[age>=30 & age<40] <- "30's"
agecat[age>=40 & age<50] <- "40's"
agecat[age>=50 & age<60] <- "50's"
agecat[age>=60 & age<70] <- "60's"
agecat[age>=70 & age<80] <- "70's"
agecat[age>=80] <- "80's"
table(agecat)

age5c<-unclass(AGE5c)	## convert into numeric
summary(age5c)
class(AGE5c)
class(age5c)

table(race)

table(sex, race)
CrossTable(sex, race)

########################################
############ Do Exercise 2B ##############
########################################

#Save your R script!

#optional, save the nhanes object you modified.
save(nhanes3_v2, file="nhanes3_20170717.rda")

## or you can save the workspace (includes all objects) with the name of your file
save.image(file="epid798.RData")

## exit R
## if you close R, you will be asked to save your workspace image



