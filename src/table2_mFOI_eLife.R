library(MASS)
library(tidyverse)
library(readstata13)
library(ggplot2)
library(gtsummary)
library(geepack)


sum_df <- read_delim("~/Desktop/Duration_sex_paper/sex_based_differences/data/table2/table2_clones_3skips.tds",  "\t")

sum_df_events <- read_delim("~/Desktop/Duration_sex_paper/sex_based_differences/data/table2/table2_events_3skips.tds",  "\t")

##factor gender correctly 
sum_df$gender = factor(sum_df$gender)
sum_df$gender = relevel(sum_df$gender, ref = "Male")

#negative binomial for overall clones, by gender
m2 <- geeglm(data = sum_df, num_clones ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male overall
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))

#IRR female to male, adjusted by age 
m6 <- geeglm(data = sum_df, num_clones ~ gender + ordered(agecat) + offset(log(days)), id = cohortid, family = "poisson")
summary(m6)

mIRR <- exp(m6$coefficients[2]) 
m_upCI <- exp(m6$coefficients[2]+1.96*(summary(m6)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m6$coefficients[2]-1.96*(summary(m6)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))


agecat1 <- sum_df %>% filter(agecat == "< 5 years")
agecat2 <- sum_df %>% filter(agecat == "5-15 years")
agecat3 <- sum_df %>% filter(agecat == "16 years or older")

#negative binomial for clones, agecat < 5 
m2 <- geeglm(data = agecat1, num_clones ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))


#negative binomial for clones, 5-15 years
m2 <- geeglm(data = agecat2, num_clones ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))


#negative binomial for clones, agecat > 16 years
m2 <- geeglm(data = agecat3, num_clones ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

##overall by agecat (no gender involved here) 
m2 <- geeglm(data = sum_df, num_clones ~ agecat + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)

exp(m2$coefficients[1])*365   #children < 5, ppy 
exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365
exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

exp(m2$coefficients[1]+m2$coefficients[2])*365   #children 5-15
exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #5-15 upper bound
exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #5-15 lower bound 

exp(m2$coefficients[1]+m2$coefficients[3])*365   #children > 16 
exp((m2$coefficients[1]+m2$coefficients[3])+1.96*(summary(m2)$coefficients[3, 2]))*365 # >16 upper bound
exp((m2$coefficients[1]+m2$coefficients[3])-1.96*(summary(m2)$coefficients[3, 2]))*365  # >16  lower bound 

#overall mFOI by clone 
m2 <- geeglm(data = sum_df, num_clones ~ num_clones + offset(log(days)), id = cohortid, family = "poisson")
summary(m2)

exp(m2$coefficients[1])*365   
exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365
exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

##############BY INFECTION EVENT###############
#FOI event
FOI_event <- labels_events %>% filter(infection_event == TRUE) %>%
  group_by(cohortid, date) %>% 
  mutate(infection_event_num = 1) %>% summarize(FOI_event = max(infection_event_num))

#join in FOI by event 
expanded <- expanded %>% left_join(FOI_event)

sum_df_events <- expanded %>% group_by(cohortid, agecat, gender) %>% 
  summarize(days = sum(days), num_events = sum(FOI_event, na.rm = TRUE)) %>% 
  mutate(mFOI = num_events/days, mFOIppy = mFOI*365)

#negative binomial for overall events, by gender
m2 <- geeglm(data = sum_df_events, num_events ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#IRR female to male overall
mIRR <- exp(m2$coefficients[2]) 
m_upCI <- exp(m2$coefficients[2]+1.96*(summary(m2)$coefficients[2, 2]))  #male upper bound CI
m_lowCI <- exp(m2$coefficients[2]-1.96*(summary(m2)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))


agecat1 <- sum_df_events %>% filter(agecat == "< 5 years")
agecat2 <- sum_df_events %>% filter(agecat == "5-15 years")
agecat3 <- sum_df_events %>% filter(agecat == "16 years or older")


#negative binomial for events, agecat < 5 
m2 <- geeglm(data = agecat1, num_events ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))


#negative binomial for events, 5-15 years
m2 <- geeglm(data = agecat2, num_events ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))

#negative binomial for events, agecat > 16 years
m2 <- geeglm(data = agecat3, num_events ~ gender + offset(log(days)), id = cohortid, family = "poisson")

mIRR <- exp(m2$coefficients[1])*365   #male incidence, ppy 
m_upCI <- exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365  #male upper bound CI
m_lowCI <- exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365
fIRR <- exp(m2$coefficients[1]+m2$coefficients[2])*365   #female incidence, ppy 
f_upCI <- exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #female upper bound
f_lowCI <- exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #female lower bound 

print(paste("Male", "inc=", round(mIRR, 2), "CI=", "(", round(m_lowCI,2), ",", round(m_upCI,2), ")", sep = " "))
print(paste("Female", "inc=", round(fIRR, 2), "CI=", "(", round(f_lowCI,2), ",", round(f_upCI,2), ")", sep = " "))


##overall by agecat (no gender involved here), negative binomial  
m2 <- geeglm(data = sum_df_events, num_events ~ agecat + offset(log(days)), id = cohortid, family = "poisson")

summary(m2)

exp(m2$coefficients[1])*365   #children < 5, ppy 
exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365
exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365

exp(m2$coefficients[1]+m2$coefficients[2])*365   #children 5-15
exp((m2$coefficients[1]+m2$coefficients[2])+1.96*(summary(m2)$coefficients[2, 2]))*365 #5-15 upper bound
exp((m2$coefficients[1]+m2$coefficients[2])-1.96*(summary(m2)$coefficients[2, 2]))*365  #5-15 lower bound 

exp(m2$coefficients[1]+m2$coefficients[3])*365   #children > 16 
exp((m2$coefficients[1]+m2$coefficients[3])+1.96*(summary(m2)$coefficients[3, 2]))*365 # >16 upper bound
exp((m2$coefficients[1]+m2$coefficients[3])-1.96*(summary(m2)$coefficients[3, 2]))*365  # >16  lower bound 

#overall mFOI by events 
m2 <- geeglm(data = sum_df_events, num_events ~ num_events + offset(log(days)), id = cohortid, family = "poisson")

summary(m2)

exp(m2$coefficients[1])*365   
exp(m2$coefficients[1]+1.96*(summary(m2)$coefficients[1, 2]))*365
exp(m2$coefficients[1]-1.96*(summary(m2)$coefficients[1, 2]))*365


#overall trend by sex (not adjusted for age), infxn event 
m4 <- geeglm(data = sum_df_events, num_events ~ gender + offset(log(days)), id = cohortid, family = "poisson")
summary(m4)$coefficients

mIRR <- exp(m4$coefficients[2])   #IRR
m_upCI <- exp(m4$coefficients[2]+1.96*(summary(m4)$coefficients[2, 2]))  #upper bound CI
m_lowCI <- exp(m4$coefficients[2]-1.96*(summary(m4)$coefficients[2, 2])) #lower bound CI 

print(paste("IRR=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))


#overall p-values, mFOI (infxn event) by sex adjusting for ordered age cat 
m6 <- geeglm(data = sum_df_events, num_events ~ gender + ordered(agecat) + offset(log(days)), id = cohortid, family = "poisson")
summary(m6)

mIRR <- exp(m6$coefficients[2])   #IRR
m_upCI <- exp(m6$coefficients[2]+1.96*(summary(m6)$coefficients[2, 2]))  #upper bound CI
m_lowCI <- exp(m6$coefficients[2]-1.96*(summary(m6)$coefficients[2, 2])) #lower bound CI 
print(paste("IRR=", round(mIRR, 3), "CI=", "(", round(m_lowCI,3), ",", round(m_upCI,3), ")", sep = " "))

#to perform sensitvity analysis, use same R code but load in different dataframes as follows:  
#2 skips
labels_clones <- read_delim("~/Desktop/ama1/src/labels_newest_clones_2skip.tab", 
                            "\t", escape_double = FALSE, col_types = cols(X1 = col_skip(), 
                                                                          burnin = col_date(format = "%Y-%m-%d"), 
                                                                          date = col_date(format = "%Y-%m-%d"), 
                                                                          enrolldate = col_date(format = "%Y-%m-%d")), 
                            trim_ws = TRUE)

labels_events <- read_delim("~/Desktop/ama1/src/labels_newest_event_2skip.tab", 
                            "\t", escape_double = FALSE, col_types = cols(X1 = col_skip(), 
                                                                          burnin = col_date(format = "%Y-%m-%d"), 
                                                                          date = col_date(format = "%Y-%m-%d"), 
                                                                          enrolldate = col_date(format = "%Y-%m-%d")), trim_ws = TRUE)



#1 skip
labels_clones <- read_delim("~/Desktop/ama1/src/labels_newest_clones_1skip.tab", 
                            "\t", escape_double = FALSE, col_types = cols(X1 = col_skip(), 
                                                                          burnin = col_date(format = "%Y-%m-%d"), 
                                                                          date = col_date(format = "%Y-%m-%d"), 
                                                                          enrolldate = col_date(format = "%Y-%m-%d")), 
                            trim_ws = TRUE)

labels_events <- read_delim("~/Desktop/ama1/src/labels_newest_event_1skip.tab", 
                            "\t", escape_double = FALSE, col_types = cols(X1 = col_skip(), 
                                                                          burnin = col_date(format = "%Y-%m-%d"), 
                                                                          date = col_date(format = "%Y-%m-%d"), 
                                                                          enrolldate = col_date(format = "%Y-%m-%d")), trim_ws = TRUE)




