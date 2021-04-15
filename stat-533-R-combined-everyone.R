rm(list=ls(all=TRUE))
options(digits=3)

#### LIBRARY USED#################
library(survival)
library(forcats)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(finalfit)
library(cmprsk)
library(MASS) 
##################################

##### LOADING DATA & SETUP#################
# if in R studio just Import Dataset->From Text don't have to use code below
if(FALSE){
  path_file <- "Melanoma.csv"
  melanoma<-read.csv(path_file,header=TRUE)
}

#getting important stats for components
ff_glimpse(melanoma)
############### SETUP ###################

#create data frame
if(FALSE) {
  "time
  Survival time in days since the operation, possibly censored.
  
  status
  The patients status at the end of the study. 1 indicates that they had died 
  from melanoma, 2 indicates that they were still alive and 3 indicates that 
  they had died from causes unrelated to their melanoma.
  
  sex
  The patients sex; 1=male, 0=female.
  
  age
  Age in years at the time of the operation.
  
  year
  Year of operation.
  
  thickness
  Tumour thickness in mm.
  
  ulcer
  Indicator of ulceration; 1=present, 0=absent."
}

melanoma.df = melanoma %>%mutate(
  
    # general motality where 2 is alive, and  1/3 are death
    status_general = ifelse(melanoma$status == 2, 0, 1), 
    
    # MM mortality comparing with 1 die by MM, 2/3 are alive
    status_mm = case_when(melanoma$status==3 ~  0, melanoma$status==2 ~ 0,melanoma$status==1 ~ 1),
    
    # Censoring adjust for competing risks 2 is alive, 1 is death by MM adjusting for 3
    status_cen = case_when(melanoma$status == 2 ~ 0, melanoma$status == 1 ~ 1, TRUE ~ 2),
    
    # labeling components for plots and tables
    thickness = ff_label(melanoma$thickness, "Tumour thickness (mm)"),
    age = ff_label(melanoma$age, "Age (Days)"), 
    sex = factor(melanoma$sex) %>% fct_recode("Male" = "1", "Female" = "0") %>% ff_label("Gender"),
    ulcer = factor(melanoma$ulcer) %>% fct_recode("No" = "0","Yes" = "1") %>% ff_label("Ulceration")
)

##################################

##### Kaplan Meier Analysis ################

# Kaplan-Meier estimates 
km_fit <- survfit(
  Surv(time/365, status_general)~1, 
  data=melanoma.df
)
summary(km_fit, times = c(1,30,60,90*(1:10)))

plot(km_fit, xlab="Days", ylab="Cummulative Probability of survival over time for MM" ,main = "Kaplan Meier plot for the survival of Maglinant Melanoma(MM)")

# BY GENDER
km_fit_gender <- survfit(
  Surv(time, status_general)~sex, 
  data=melanoma.df
)

xlabel = "Days"
ylabel = "Cummulative Probability of survival over time for MM"
plot_title = "Kaplan Meier plot for the survival of Maglinant Melanoma(MM) by gender"
km_gender_plot <- autoplot(km_fit_gender)
km_gender_plot <- km_gender_plot + 
  ggtitle(plot_title) +
  labs(x = xlabel, y = ylabel) +
  guides(fill=FALSE) +
  labs(colour = "Gender") +
  scale_color_manual(labels = c("Female", "Male"), values = c(1, 2))
print(km_gender_plot)

# BY ULCER
km_fit_ulcer <- survfit(
  Surv(time, status_general)~ulcer, 
  data=melanoma.df
)

plot_title = "Kaplan Meier plot for the survival of Maglinant Melanoma(MM) by ulceration"
km_ulcer_plot <- autoplot(km_fit_ulcer)
km_ulcer_plot <- km_ulcer_plot + 
  ggtitle(plot_title) +
  labs(x = xlabel, y = ylabel) +
  guides(fill=FALSE) +
  labs(colour = "Ulceration") +
  scale_color_manual(labels = c( "Ulcer Absent", "Ulcer Present"), values = c(1, 2))
print(km_ulcer_plot)


# BY AGES

# Let cutt off age be 50 years of age
# from https://www.munichre.com/content/dam/munichre/marc/pdf/malignant-melanoma-may-2020.pdf/_jcr_content/renditions/original./malignant-melanoma-may-2020.pdf
if(FALSE){
  "In the U.S., the lifetime risk of developing melanoma was estimated to be about 1 in 28 for
  Caucasians in the year 2020. The incidence of melanoma is relatively low in dark-skinned
  populations. Overall, the lifetime risk of developing melanoma is greater for men than
  women; however, for those under age 50, the inverse is true. The incidence rate begins to
  rise after puberty, increases until the age of 65 to 70 years, and then decreases with the
  average age at diagnosis being around 65." 
}

plot_km_age <- function(cutoff_age_arr) {
  for (a in cutoff_age_arr) {
    lt_age = paste("Less than ", a, sep="")
    ov_age = paste("Over ", a, sep="")
    mm_age_adj = melanoma.df %>% mutate(AG = ifelse((age < a), lt_age, ov_age),
                 AG = factor(AG))
    
    km_fit_age <- survfit(
      Surv(time, status_general)~AG, 
      data=mm_age_adj
    )
    
    plot_title = "Kaplan Meier plot for the survival of Maglinant Melanoma(MM) by age"
    km_age_plot <- autoplot(km_fit_age)
    km_age_plot <- km_age_plot + 
      ggtitle(plot_title) +
      labs(x = xlabel, y = ylabel) +
      guides(fill=FALSE) +
      labs(colour = "AGE") +
      scale_color_manual(labels = c(lt_age, ov_age), values = c(1, 2))
    print(km_age_plot) 
  }
}

cutoff_ages <- c(40,50,60,70)
plot_km_age(cutoff_ages)
####################################################################################
######## Cox Prop Haz Model
#model fit and variable selection 
# Perform CoxPropHaz model
modelfit.cox.ph = coxph(Surv(time,death_status)~strata(factor(ulcer))+thickness+age+factor(sex)+year, method=c("breslow"), data=melanoma.df)
modelfit.cox.ph
summary(modelfit.cox.ph)

# testing proportional hazard assumption

#gerneral
par(mfrow=c(2,3))
death_status <-  melanoma.df$status_general
fit.cox = coxph(Surv(time,death_status)~factor(sex)+age+thickness+year+factor(ulcer),method=c("breslow"), na.action=na.exclude, data=melanoma.df)
fit.cox.zph = cox.zph(fit.cox)
summary(fit.cox.zph)
plot(fit.cox.zph)

fit.cox.zph
# A correlation of zero indicates that the model met the proportional hazards assumption
#             chisq df     p
#factor(sex)   0.503  1 0.478
#age           2.063  1 0.151
#thickness     2.825  1 0.093 
#year          0.450  1 0.502
#factor(ulcer) 4.320  1 0.038 --- reject Ho, PH does not hold
#GLOBAL        7.870  5 0.164

# we need to stratify ulceration component
par(mfrow=c(2,3))
fit.cox = coxph(Surv(time,death_status)~factor(sex)+age+thickness+year+strata(factor(ulcer)),method=c("breslow"), na.action=na.exclude, data=melanoma.df)
fit.cox.zph = cox.zph(fit.cox)
summary(fit.cox.zph)
plot(fit.cox.zph)
abline(h=0, lty=3)
fit.cox.zph

#model testing and model selection

#perform local tests for each explanatory factor in the model
#We know thickness and ulcer factor are absolute effective but we need to once test age, gender and year 
#to figure out more effective factors.
#will use Wald test 
coxph.localtest<-function(coxph.fit, df){
  coef<-coxph.fit$coef
  var<-coxph.fit$var
  loglik<-coxph.fit$loglik[2]
  p<-length(coef)
  AIC<- -2*loglik+2*p    #Using the formula on p.277;
  var.mat<-solve(as.matrix(var[(p-df+1):p, (p-df+1):p]))
  coe.mat<-matrix(c(coef[(p-df+1):p]), df, 1)
  WaldChiSq<-t(coe.mat)%*%var.mat%*%coe.mat
  pvalue<-1-pchisq(WaldChiSq,df)
  results<-c(df, WaldChiSq, pvalue, AIC)
  list(results)
}

options(width=60,length=200,digits=5) 

Table_checkup<-matrix(0,3,4)

coxph.sex<-coxph(Surv(time,death_status)~strata(factor(ulcer))+thickness+factor(sex), method=c("breslow"), data=melanoma.df)
coxph.age<-coxph(Surv(time,death_status)~strata(factor(ulcer))+thickness+age, method=c("breslow"), data=melanoma.df)
coxph.year <- coxph(Surv(time,death_status)~strata(factor(ulcer))+thickness+year, method=c("breslow"), data=melanoma.df)

Table_checkup[1,]<-c(coxph.localtest(coxph.sex, df=1)[[1]])
Table_checkup[2,]<-c(coxph.localtest(coxph.age, df=1)[[1]])
Table_checkup[3,]<-c(coxph.localtest(coxph.year, df=1)[[1]])

cat("Table 8.6: Local test for possible confounders, adjusted for factors of gender, age, thickness and ulceration", "\n")
cat("    DF, Wald Chi-Squre, p-value, AIC:", "\n")
print(Table_checkup)

#since all p_value except age more than 0.05, we fail to reject Ho and said sex and year is not effective as only age is 

#This is the finalized model
coxph_final<-coxph(Surv(time,death_status)~strata(factor(ulcer))+thickness+age, method=c("breslow"), data=melanoma.df)

cat("Table 8.9: Analysis of Variance Table for the Final Model for the survival of Malignant Melanoma", "\n")
print(coxph_final)

######################################
#Double check it by taking forward and backward AIC method but it turn out same as reject age, year and gender

fit.cox <- coxph(Surv(time,death_status)~1,data=melanoma.df)
fit.cox

# backward AIC method
stepAIC(fit.cox,  ditrection="backward",trace=-1) 

# forward AIC method
stepAIC(fit.cox, direction="forward", 
        scope=list(lower=fit.cox, upper=~factor(sex)+age+thickness+year+strata(factor(ulcer))))
cat("with backward and forward AIC method, we obtain the same results as hazard proportion.")
#######################################
# Adjust for competing risks 

explanatory = c("age", "sex", "thickness", "ulcer")
mm_class = "Surv(time, status_mm)"
cen_class = "Surv(time, status_cen)"

melanoma.df %>%coxphmulti(mm_class, explanatory)
melanoma.df%>%crrmulti(cen_class, explanatory)
##################################

cat("Due to our cox hazard proportion and Kaplan Meier Analysis, we can conlude that thickness, ulcer status and age are
    effective to death rate or status, since only male is effective but not female, then we also conclude gender and 
    operation year is not effective factors.")
