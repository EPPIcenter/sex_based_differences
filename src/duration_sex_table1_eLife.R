library(MASS)
library(tidyverse)
library(ggplot2)
library(EnvStats)
library(geepack)

#to get person-time 

#load in incidence databases 
overall_inc <- read_delim("~/Desktop/Duration_sex_paper/sex_based_differences/data/table1/overall_incidence_table1.tds", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) 

agecat_inc <- read_delim("~/Desktop/Duration_sex_paper/sex_based_differences/data/table1/agecat_inc.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) 

exp1 <- read_delim("~/Desktop/Duration_sex_paper/sex_based_differences/data/table1/llin_dates.tsv", 
                   "\t", escape_double = FALSE, trim_ws = TRUE) 


#load in full_meta databases
allvisits <- read_delim("~/Desktop/ama1/prism2/stata/full_meta_6mo.tab", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)


overall_inc$gender = factor(overall_inc$gender)
overall_inc$gender = relevel(overall_inc$gender, ref = "Male")

agecat_inc$gender = factor(agecat_inc$gender)
agecat_inc$gender = relevel(agecat_inc$gender, ref = "Male")

allvisits$gender = factor(allvisits$gender)
allvisits$gender = relevel(allvisits$gender, ref = "Male")

#number of people
age_sex_person <- allvisits %>% arrange(date) %>% group_by(cohortid) %>% 
  slice(1) %>% select(cohortid, gender, agecat) %>% unique()

sex <- age_sex_person %>% select(cohortid, gender)

age_sex_person %>% filter(gender == "Male") #233
age_sex_person %>% filter(gender == "Female") #244

age_sex_person %>% filter(agecat == "< 5 years" & gender == "Male")  #73

age_sex_person %>% filter(agecat == "< 5 years" & gender == "Female") #84

age_sex_person %>% filter(agecat == "5-15 years" & gender == "Male")  #101
age_sex_person %>% filter(agecat == "5-15 years" & gender == "Female") #71

age_sex_person %>% filter(agecat == "16 years or older" & gender == "Male")  #59
age_sex_person %>% filter(agecat == "16 years or older" & gender == "Female") #89

#median time 
time_per_person <- agecat_inc %>% group_by(cohortid, gender) %>% summarize(days = sum(days))

time_per_person %>% filter(gender == "Male") %>% summary()
time_per_person %>% filter(gender == "Female") %>% summary()

overall_inc <- overall_inc %>% left_join(age_sex_person)

overall_inc %>% filter(agecat == "< 5 years" & gender == "Male") %>% summary()
overall_inc %>% filter(agecat == "< 5 years" & gender == "Female") %>% summary()

overall_inc %>% filter(agecat == "5-15 years" & gender == "Male")  %>% summary()
overall_inc %>% filter((agecat == "5-15 years" & gender == "Female") ) %>% summary()

overall_inc %>% filter(agecat == "16 years or older" & gender == "Male") %>% summary()
overall_inc %>% filter(agecat == "16 years or older" & gender == "Female") %>% summary()

#LLIN use 

#overall by gender
exp1 %>% filter(gender == "Male") %>% summarize(percent = sum(llin)/sum(days)) 
exp1 %>% filter(gender == "Female") %>% summarize(percent = sum(llin)/sum(days)) 

#children < 5 
exp1 %>% filter(agecat == "< 5 years" & gender == "Male") %>% summarize(percent = sum(llin)/sum(days)) 
exp1 %>% filter(agecat == "< 5 years" & gender == "Female") %>% summarize(percent = sum(llin)/sum(days)) 

#children 5-15 
exp1 %>% filter(agecat == "5-15 years" & gender == "Male") %>% summarize(percent = sum(llin)/sum(days)) 
exp1 %>% filter(agecat == "5-15 years" & gender == "Female") %>% summarize(percent = sum(llin)/sum(days)) 

#adults (16+)
exp1 %>% filter(agecat == "16 years or older" & gender == "Male") %>% summarize(percent = sum(llin)/sum(days)) 
exp1 %>% filter(agecat == "16 years or older" & gender == "Female") %>% summarize(percent = sum(llin)/sum(days)) 

#malaria incidence by sex 
m2 <- geeglm(data = overall_inc, anymalaria ~ gender + offset(log(days)), id = cohortid, family = "poisson")
mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male, not adjusted by age
#make sure reference is male
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))

##IRR female to male for malaria by sex, adjusted for age 
y <- geeglm(data = agecat_inc, anymalaria ~ gender + agecat + offset(log(days)), id = cohortid, family = "poisson")
summary(y)

mIRR <- exp(y$coefficients[2])   #IRR
m_upCI <- exp(y$coefficients[2]+1.96*(summary(y)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(y$coefficients[2]-1.96*(summary(y)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

#malaria incidence by sex/agecat 
#<5 years old
expanded_agecat1 <- agecat_inc %>% filter(agecat == "< 5 years")
m2 <- geeglm(data = expanded_agecat1, anymalaria ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#5-15 years
expanded_agecat2 <- agecat_inc %>% filter(agecat == "5-15 years")
m2 <- geeglm(data = expanded_agecat2, anymalaria ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#16 years and greater 
expanded_agecat3 <- agecat_inc %>% filter(agecat == "16 years or older")
m2 <- geeglm(data = expanded_agecat3, anymalaria ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#TRAVEL 

#travel incidence by sex 
m2 <- geeglm(data = overall_inc, travelout ~ gender + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)
mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#travel incidence by sex/agecat 
#<5 years old
m2 <- geeglm(data = expanded_agecat1, travelout ~ gender + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)
mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male, agecat1
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))

#5-15 years
m2 <- geeglm(data = expanded_agecat2, travelout ~ gender + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)
mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male, agecat2
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))

#16 years and greater 

m2 <- geeglm(data = expanded_agecat3, travelout ~ gender + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)
mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male, agecat3 (>16)
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))


#Now by visit
allvisits1 <- allvisits %>% 
  filter(visittype == "routine visit" | visittype == "enrollment visit") %>% 
  filter(date <= "2019-03-31")

allvisits1 <- allvisits1 %>% 
  mutate(qpcr_one = if_else(qpcr > 0, 1, 0), LM_one=ifelse(parasitedensity > 0, 1, 0)) %>% 
  mutate(itn = ifelse(itnlastnight == "no", 1, 0)) %>% 
  mutate(travel_one = ifelse(travelout == "yes", 1, 0)) %>% 
  mutate(qpcr_one = ifelse(is.na(qpcr_one), 0, qpcr_one)) %>% 
  mutate(LM_one = ifelse(is.na(LM_one), 0, LM_one)) %>% 
  mutate(travel_one = ifelse(is.na(travel_one), 0, travel_one)) %>% 
  mutate(n = 1)

#how many routine/enrollment visits in total
#by sex & agecat
visits_total <- allvisits1 %>% 
  group_by(agecat, gender) %>% summarize(visits = sum(n()))

#by sex
visits_total1 <- allvisits1 %>% 
  group_by(gender) %>% summarize(visits = sum(n()))

##QPCR prevalence  
#significant, p = 0.03  
x<- geeglm(qpcr_one ~ gender, data = allvisits1, id = cohortid, family = "poisson")
summary(x)

mIRR <- exp(x$coefficients[1])   #male proportion, ppy 
m_upCI <- exp(x$coefficients[1]+1.96*(summary(x)$coefficients[1, 2]))  #male upper bound CI

m_lowCI <- exp(x$coefficients[1]-1.96*(summary(x)$coefficients[1, 2])) #lower bound CI 
fIRR <- exp(x$coefficients[1]+x$coefficients[2])   #female proportion, ppy 
f_upCI <- exp((x$coefficients[1]+x$coefficients[2])+1.96*(summary(x)$coefficients[2, 2])) #female upper bound
f_lowCI <- exp((x$coefficients[1]+x$coefficients[2])-1.96*(summary(x)$coefficients[2, 2]))  #female lower bound

print(paste("Male", "prop=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))
print(paste("Female", "prop=", round(fIRR, 3), "CI=", "(", round(f_lowCI,3), ",", round(f_upCI,3), ")", sep = " "))

#prevalence ratio female to male, not adjusted by age 
mIRR <- exp(summary(x)$coefficients[2,1]) 
m_upCI <- exp(summary(x)$coefficients[2,1]+1.96*(summary(x)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(summary(x)$coefficients[2,1]-1.96*(summary(x)$coefficients[2, 2])) #lower bound CI 

print(paste("prev ratio=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

#qPCR less than 5 years 
x <- geeglm(qpcr_one ~ gender, data = allvisits1 %>% filter(agecat == "< 5 years"), id = cohortid, family = poisson)
summary(x)

#qpcr 5-15
x <- geeglm(qpcr_one ~ gender, data = allvisits1 %>% filter(agecat == "5-15 years"), id = cohortid, family = poisson)
summary(x)

#qpcr >=16
x <- geeglm(qpcr_one ~ gender, data = allvisits1 %>% filter(agecat == "16 years or older"), id = cohortid, family = poisson)
summary(x)

####light microscopy, overall 
x <- geeglm(LM_one ~ gender, data = allvisits1, id = cohortid, family = poisson)
summary(x)

mIRR <- exp(x$coefficients[1])   #male proportion, ppy 
m_upCI <- exp(x$coefficients[1]+1.96*(summary(x)$coefficients[1, 2]))  #male upper bound CI
m_lowCI <- exp(x$coefficients[1]-1.96*(summary(x)$coefficients[1, 2])) #lower bound CI 
fIRR <- exp(x$coefficients[1]+x$coefficients[2])   #female proportion, ppy 
f_upCI <- exp((x$coefficients[1]+x$coefficients[2])+1.96*(summary(x)$coefficients[2, 2])) #female upper bound
f_lowCI <- exp((x$coefficients[1]+x$coefficients[2])-1.96*(summary(x)$coefficients[2, 2]))  #female lower bound

print(paste("Male", "prop=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))
print(paste("Female", "prop=", round(fIRR, 3), "CI=", "(", round(f_lowCI,3), ",", round(f_upCI,3), ")", sep = " "))

#prevalence ratio female to male, not adjusted by age 
mIRR <- exp(summary(x)$coefficients[2,1]) 
m_upCI <- exp(summary(x)$coefficients[2,1]+1.96*(summary(x)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(summary(x)$coefficients[2,1]-1.96*(summary(x)$coefficients[2, 2])) #lower bound CI 

print(paste("prev ratio=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

#LM <5
x <- geeglm(LM_one ~ gender, data = allvisits1 %>% filter(agecat == "< 5 years"), id = cohortid, family = poisson)
summary(x)

#LM 5-15
x <- geeglm(LM_one ~ gender, data = allvisits1 %>% filter(agecat == "5-15 years"), id = cohortid, family = poisson)
summary(x)

#LM > 16 years
x <- geeglm(LM_one ~ gender, data = allvisits1 %>% filter(agecat == "16 years or older"), id = cohortid, family = poisson)
summary(x)


#qpcr density by gender, adjusting for age
#all visits including malaria visits
x <- geeglm( log10(qpcr) ~ gender + ageyrs, data = allvisits %>% filter(qpcr >0), id = cohortid)
summary(x)

#all routine/enrollment visits (asymptomatic)
x <- geeglm( log10(qpcr) ~ gender + ageyrs, data = allvisits1 %>% filter(qpcr >0), id = cohortid)
summary(x)

#geometric mean parasite density by sex of routine/enrollment visits 

allvisits1 %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Male") %>% 
  summarize(geoMean(qpcr))


allvisits1 %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Female") %>% 
  summarize(geoMean(qpcr))

#geometric mean parasite density, agecat1 

agecat1_routines <- allvisits1 %>% filter(agecat == "< 5 years")

agecat1_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Male") %>% 
  summarize(geoMean(qpcr))

agecat1_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Female") %>% 
  summarize(geoMean(qpcr))

#geometric mean parasite density, agecat2

agecat2_routines <- allvisits1 %>% filter(agecat == "5-15 years")

agecat2_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Male") %>% 
  summarize(geoMean(qpcr))

agecat2_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Female") %>% 
  summarize(geoMean(qpcr))

#geometric mean parasite density, agecat3

agecat3_routines <- allvisits1 %>% filter(agecat == "16 years or older")

agecat3_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Male") %>% 
  summarize(geoMean(qpcr))

agecat3_routines %>% 
  select(cohortid, date, agecat, gender, qpcr) %>% 
  filter(qpcr > 0 & !is.na(qpcr) & gender == "Female") %>% 
  summarize(geoMean(qpcr))


#MOI 
interval_clone <- read_tsv("~/Desktop/Duration_sex_paper/sex_based_differences/data/labeled_db/newsest_labels_clones.tsv")
MOI <- interval_clone %>% select(cohortid, date, h_popUID) %>% group_by(cohortid, date) %>%
  summarize(MOI = n())

interval_clone <- interval_clone %>% left_join(MOI)

#MOI overall
interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Male") %>% 
  summary()

interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Female") %>% 
  summary()


#MOI age 5 or less 
interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Male" & agecat == "< 5 years") %>% 
  summary()

interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Female" & agecat == "< 5 years") %>% 
  summary()

#MOI 5-15
interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Male" & agecat == "5-15 years") %>% 
  summary()

x <- interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Male" & agecat == "5-15 years")

hist(x$MOI)

interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Female" & agecat == "5-15 years") %>% 
  summary()

y <- interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Female" & agecat == "5-15 years")

hist(y$MOI)

#MOI 16+
interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Male" & agecat == "16 years or older") %>% 
  summary()

interval_clone %>% 
  select(cohortid,agecat, gender, MOI) %>% 
  filter(gender == "Female" & agecat == "16 years or older") %>% 
  summary()


###########Statistical test showing sex is a significant predictor of prevalence 

table(allvisits1$osantimalarial)

#can't use osantimalarial as a predictor, won't converge 

exp2 <- exp1 %>% 
  group_by(cohortid, monthyear) %>% summarize(llin_total = sum(llin, na.rm = TRUE)/sum(days))

allvisits1 <- allvisits1 %>% left_join(exp2) 
allvisits1 <- allvisits1 %>% mutate(llin_total = ifelse(is.na(llin_total), 1, llin_total)) 

m3 <- geeglm(data = allvisits1, qpcr_one ~ llin_total + travel_one + agecat + gender, id = cohortid, family = "poisson" )
summary(m3)

#prevalence ratio female to male
mIRR <- exp(summary(m3)$coefficients[6,1]) 
m_upCI <- exp(summary(m3)$coefficients[6,1]+1.96*(summary(m3)$coefficients[6, 2]))  # upper bound CI
m_lowCI <- exp(summary(m3)$coefficients[6,1]-1.96*(summary(m3)$coefficients[6, 2])) #lower bound CI 

print(paste("prev ratio=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

m3 <- geeglm(data = allvisits1, LM_one ~ travel_one + llin_total + agecat + gender, id = cohortid, family = "poisson" )
summary(m3)

#prevalence ratio female to male, not adjusted by age + llin + travel 
mIRR <- exp(summary(m3)$coefficients[6,1]) 
m_upCI <- exp(summary(m3)$coefficients[6,1]+1.96*(summary(m3)$coefficients[6, 2]))  # upper bound CI
m_lowCI <- exp(summary(m3)$coefficients[6,1]-1.96*(summary(m3)$coefficients[6, 2])) #lower bound CI 

print(paste("prev ratio=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

