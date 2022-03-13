##########################################
##########################################
############# Normal LASSO ###############
##########################################
##########################################
library(tidyverse)
library(MatchIt)
library(twang)
library(glmnet)
library(ggplot2)
library(dplyr)
library(lsr)
library(survey)
library(rbounds)
library(naniar)
library(tableone)
library(optmatch)
library(mice)
library(plotmo)
library(corrplot)
library(ggplot2)
library(forestplot)
library(ggcorrplot)
library(MatchThem)
library(gtsummary)
library(rstatix)
setwd("/Users/Crystal/Desktop/Practicum/Data")
NAFLD <-read_csv("data_NAFLD_PS.csv",show_col_types = FALSE)
NAFLD<-NAFLD %>%mutate(Sex=as.factor(ifelse(Sex=="Male", 0, 1)))
NAFLD<-NAFLD[,c(3:4,7,8,10:12,14,27,41)]
variables<-c("Sex","LDL","Triglycerides","Glycated_haemoglobin__HbA1c_",
             "HDL_cholesterol","Waist_circumference","Systolic_Blood_Pressure","BMIdich")

################################
################################
############# CCA ###############
################################
################################

omitNAFLD<- na.omit(NAFLD)
newvariables<-c("LDL","Triglycerides","Glycated_haemoglobin__HbA1c_","Waist_circumference","Systolic_Blood_Pressure","BMIdich")
ps.formula <- as.formula((paste("T2D ~",paste(newvariables,collapse = "+"))))
param<-matchit(ps.formula,data = omitNAFLD)
omitNAFLD$param_ps<-param$distance
set.seed(123456)
nonparam<-ps(as.formula((paste("T2D ~ ",paste(newvariables,collapse = "+")))),
             data = as.data.frame(omitNAFLD),n.trees=1000,interaction.depth=9,shrinkage=0.01,
             stop.method="es.mean",estimand="ATT")
omitNAFLD$nonparam_ps<-nonparam$ps$es.mean
summary(nonparam$gbm.obj)
ggplot(data = omitNAFLD,aes(param_ps,fill=factor(T2D)))+
  geom_histogram(binwidth = 0.005,alpha = 0.5,position="identity")
logitPS<-sd(-log(1/omitNAFLD$param_ps-1))

## Greedy matching

greedy<-matchit(ps.formula,data=omitNAFLD,
                distance="logit",
                caliper=0.1*sd(logitPS),
                method="nearest")
greedy
greedy$nn <- c(3695, 1382)


## Optimal matching
optimal <- matchit(ps.formula,data= omitNAFLD,
                   distance="logit",
                   method = "optimal")
optimal
optimal$nn <- c(3695, 1382)


## 1:2 matching
one_two <- matchit(ps.formula,data=omitNAFLD,
                   distance="logit",
                   caliper=0.1*sd(logitPS),
                   method = "nearest", ratio = 2)
one_two
one_two$nn <-c(3695, 2073)

## Full matching
full <- matchit(ps.formula,data=omitNAFLD,
                distance="logit",
                caliper=0.1*sd(logitPS),
                method = "full")
full
full$nn<- c(3695, 3695)
matching.t<-matrix(c(greedy$nn,optimal$nn,one_two$nn,full$nn), ncol=2, byrow=TRUE)
rownames(matching.t)<-c("Greedy matching","Optimal matching","1:2 matching","Full matching")
colnames(matching.t)<-c("Original Groups","Matched Groups")
knitr::kable(as.table(matching.t))

## Subclassification
omitNAFLD <- omitNAFLD %>% mutate(subclass = ntile(param_ps, 5))


## Subclassification with full matching
subc <- matchit(ps.formula, data = omitNAFLD,
                method = "subclass",
                subclass = 5)
subc
subc.t<-matrix(c(3695, 3695), ncol=2, byrow=TRUE)
colnames(subc.t)<-c("Original Groups","Matched Groups")
rownames(subc.t)<-"Subclasses"
knitr::kable(as.table(subc.t))
plot(subc, type = "jitter")

### Weighting

## IPTW with stabilization
omitNAFLD$iptw<- ifelse(omitNAFLD$T2D == 1, (mean(omitNAFLD$param_ps))/omitNAFLD$param_ps, (mean(1-omitNAFLD$param_ps))/(1-omitNAFLD$param_ps))

# NAFLD_status~T2D
formula1<-as.formula((paste("NAFLD_status","~ T2D+ ",paste(newvariables,collapse = "+"))))
## Greedy matching
matched_greedy <- match.data(greedy,subclass = "block")
glm_NAFLD_greedy <- glm(formula1,
                        data = matched_greedy,
                        family = binomial)
summary_greedy <- summary(glm_NAFLD_greedy)
exp(coef(glm_NAFLD_greedy))


## Optimal matching
matched_optimal <- match.data(optimal,subclass = "block")
glm_NAFLD_optimal <- glm(formula1,
                         data = matched_optimal,
                         family = binomial)
summary_optimal <- summary(glm_NAFLD_optimal)
exp(coef(glm_NAFLD_optimal))


## 1:2 matching
matched_one_two <- match.data(one_two,subclass = "block")
glm_NAFLD_one_two <- glm(formula1, 
                         data = matched_one_two,
                         family = binomial)
summary_one_two <- summary(glm_NAFLD_one_two)
exp(coef(glm_NAFLD_one_two))


## Full matching
matched_full <- match.data(full,subclass = "block")
glm_NAFLD_full <- glm(formula1,
                      data = matched_full,
                      family = binomial, weights = weights)
summary_full <-summary(glm_NAFLD_full)
exp(coef(glm_NAFLD_full))


## Subclassification with full matching
subced <- match.data(subc,subclass = "block")
design <- svydesign(id = ~1, strata = ~subclass, weights = ~weights, nested = TRUE, data = subced)
svyglm_NAFLD_subc <- svyglm(formula1, design = design, family = binomial)
summary_subc <- summary(svyglm_NAFLD_subc)
exp(coef(svyglm_NAFLD_subc))


## Weighting by IPTW with stabilization
glm_NAFLD_iptw <- glm(formula1,
                      data = omitNAFLD,
                      weights = iptw,
                      family = binomial)
summary_iptw <-summary(glm_NAFLD_iptw)
exp(coef(glm_NAFLD_iptw))

NAFLD.T2D.HR.est<-c()
NAFLD.T2D.logHR.est<-c()
NAFLD.T2D.se<-c()
NAFLD.T2D.pvalue<-c()
NAFLD.T2D.lowerCI<-c()
NAFLD.T2D.upperCI<-c()

NAFLD.T2D.HR.est[1] <-exp(summary_greedy$coefficients[2,1])
NAFLD.T2D.logHR.est[1]<-summary_greedy$coefficients[2,1]
NAFLD.T2D.se[1]<-summary_greedy$coefficients[2,2]
NAFLD.T2D.pvalue[1]<-summary_greedy$coefficients[2,4]
NAFLD.T2D.lowerCI[1]<- exp(summary_greedy$coefficients[2,1]-1.96*summary_greedy$coefficients[2,2])
NAFLD.T2D.upperCI[1]<- exp(summary_greedy$coefficients[2,1]+1.96*summary_greedy$coefficients[2,2])

NAFLD.T2D.HR.est[2] <-exp(summary_optimal$coefficients[2,1])
NAFLD.T2D.logHR.est[2]<-summary_optimal$coefficients[2,1]
NAFLD.T2D.se[2]<-summary_optimal$coefficients[2,2]
NAFLD.T2D.pvalue[2]<-summary_optimal$coefficients[2,4]
NAFLD.T2D.lowerCI[2]<- exp(summary_optimal$coefficients[2,1]-1.96*summary_optimal$coefficients[2,2])
NAFLD.T2D.upperCI[2]<- exp(summary_optimal$coefficients[2,1]+1.96*summary_optimal$coefficients[2,2])

NAFLD.T2D.HR.est[3] <-exp(summary_one_two$coefficients[2,1])
NAFLD.T2D.logHR.est[3]<-summary_one_two$coefficients[2,1]
NAFLD.T2D.se[3]<-summary_one_two$coefficients[2,2]
NAFLD.T2D.pvalue[3]<-summary_one_two$coefficients[2,4]
NAFLD.T2D.lowerCI[3]<- exp(summary_one_two$coefficients[2,1]-1.96*summary_one_two$coefficients[2,2])
NAFLD.T2D.upperCI[3]<- exp(summary_one_two$coefficients[2,1]+1.96*summary_one_two$coefficients[2,2])

NAFLD.T2D.HR.est[4] <-exp(summary_full$coefficients[2,1])
NAFLD.T2D.logHR.est[4]<-summary_full$coefficients[2,1]
NAFLD.T2D.se[4]<-summary_full$coefficients[2,2]
NAFLD.T2D.pvalue[4]<-summary_full$coefficients[2,4]
NAFLD.T2D.lowerCI[4]<- exp(summary_full$coefficients[2,1]-1.96*summary_full$coefficients[2,2])
NAFLD.T2D.upperCI[4]<- exp(summary_full$coefficients[2,1]+1.96*summary_full$coefficients[2,2])

NAFLD.T2D.HR.est[5] <-exp(summary_subc$coefficients[2,1])
NAFLD.T2D.logHR.est[5]<-summary_subc$coefficients[2,1]
NAFLD.T2D.se[5]<-summary_subc$coefficients[2,2]
NAFLD.T2D.pvalue[5]<-summary_subc$coefficients[2,4]
NAFLD.T2D.lowerCI[5]<- exp(summary_subc$coefficients[2,1]-1.96*summary_subc$coefficients[2,2])
NAFLD.T2D.upperCI[5]<- exp(summary_subc$coefficients[2,1]+1.96*summary_subc$coefficients[2,2])

NAFLD.T2D.HR.est[6] <-exp(summary_iptw$coefficients[2,1])
NAFLD.T2D.logHR.est[6]<-summary_iptw$coefficients[2,1]
NAFLD.T2D.se[6]<-summary_iptw$coefficients[2,2]
NAFLD.T2D.pvalue[6]<-summary_iptw$coefficients[2,4]
NAFLD.T2D.lowerCI[6]<- exp(summary_iptw$coefficients[2,1]-1.96*summary_iptw$coefficients[2,2])
NAFLD.T2D.upperCI[6]<- exp(summary_iptw$coefficients[2,1]+1.96*summary_iptw$coefficients[2,2])

NAFLDoutcome.T2D.result <- 
  structure(list(
    mean  = c(NA, NA, NAFLD.T2D.HR.est), 
    lower = c(NA, NA, NAFLD.T2D.lowerCI),
    upper = c(NA, NA, NAFLD.T2D.upperCI)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

A <- cbind2(c("","Method","Greedy Matching","Optimal Matching","1:2 Matching","Full Matching","Subclassification", "IPTW"),
            c("","Effect", as.character(round(NAFLD.T2D.HR.est,3))))
B<- cbind2(c("","SE", as.character(round(NAFLD.T2D.se,3))),
           c("","p-value", ifelse(NAFLD.T2D.pvalue<0.001, "<0.001", round(NAFLD.T2D.pvalue,3))))
tabletext <- cbind2(A ,B)

forestplot(tabletext, 
           NAFLDoutcome.T2D.result,new_page = TRUE,
           col = fpColors(box="royalblue",line="darkblue"),
           xlab = "Estimate with 95% CI",
           hrzl_lines= gpar(col = "#444444"),
           is.summary = c(TRUE, TRUE, rep(FALSE, 10)),
           clip = c(0.5,5.5), zero=0.5,
           txt_gp=fpTxtGp(label = gpar(cex = 0.8),
                          title = gpar(cex = 1),
                          ticks = gpar(cex = 0.7),
                          xlab = gpar(cex = 0.7)),
           grid = structure(seq(0.5, 5.5, by=0.25), 
                            gp = gpar(lty = 2, col = "#CCCCFF")))
grid.text("Exposure of T2D on NAFLD", .5, 0.98, gp=gpar(fontface="bold"))

################################
################################
############# MI ###############
################################
################################
#Impute the data with 5 imputations (more is better)
mi.data<-mice(NAFLD, m = 5)
c.data<-complete(mi.data)

################################
############MIte################
################################

## Greedy matching
greedy.matched.MIte.datasets <- matchthem(ps.formula, mi.data,
                                          approach = 'within',
                                          method = 'nearest',
                                          caliper = 0.05)
## Optimal matching
optimal.matched.MIte.datasets <- matchthem(ps.formula, mi.data,
                                           approach = 'within',
                                           method = 'optimal',
                                           caliper = 0.05)
## 1:2 matching
one_two.matched.MIte.datasets <- matchthem(ps.formula, mi.data,
                                           approach = 'within',
                                           method = 'nearest',
                                           caliper = 0.05, ratio = 2)
## Full matching
full.matched.MIte.datasets <- matchthem(ps.formula,mi.data,
                                        approach = 'within',
                                        caliper = 0.05,
                                        method = "full")

## Subclassification with full matching
subc.MIte.datasets <- matchthem(ps.formula, mi.data,
                                approach = 'within',
                                method = "subclass",
                                subclass = 5)
subc.MIte.datasets
plot(subc.MIte.datasets, type = "jitter")

# IPTW
weighted.MIte.datasets <- weightthem(ps.formula, mi.data,
                                     approach = 'within', 
                                     method = 'ps', estimand = 'ATM')
weighted.MIte.datasets

library(cobalt)
bal.tab(greedy.matched.MIte.datasets, abs = TRUE)
bal.tab(optimal.matched.MIte.datasets, abs = TRUE)
bal.tab(one_two.matched.MIte.datasets, abs = TRUE)
bal.tab(full.matched.MIte.datasets, abs = TRUE)
bal.tab(subc.MIte.datasets, abs = TRUE)
bal.tab(weighted.MIte.datasets, abs = TRUE)

library(survey)
formula1<-as.formula((paste("NAFLD_status","~ T2D+ ",paste(newvariables,collapse = "+"))))

greedy.matched.MIte.models <- with(data = greedy.matched.MIte.datasets,
                                   expr = svyglm(formula1, family = quasibinomial))
greedy.matched.MIte.results <- pool(greedy.matched.MIte.models)
summary_MIte_greedy <- summary(greedy.matched.MIte.results, conf.int = TRUE)

optimal.matched.MIte.models <- with(data = optimal.matched.MIte.datasets,
                                    expr = svyglm(formula1, family =quasibinomial))
optimal.matched.MIte.results <- pool(optimal.matched.MIte.models)
summary_MIte_optimal <-summary(optimal.matched.MIte.results, conf.int = TRUE)

one_two.matched.MIte.models <- with(data = one_two.matched.MIte.datasets,
                                    expr = svyglm(formula1, family = quasibinomial))
one_two.matched.MIte.results <- pool(one_two.matched.MIte.models)
summary_MIte_one_two <- summary(one_two.matched.MIte.results, conf.int = TRUE)

full.matched.MIte.models <- with(data = full.matched.MIte.datasets,
                                 expr = svyglm(formula1, family = quasibinomial))
full.matched.MIte.results <- pool(full.matched.MIte.models)
summary_MIte_full <- summary(full.matched.MIte.results, conf.int = TRUE)


subc.MIte.models <- with(data = subc.MIte.datasets,
                         expr = svyglm(formula1, family = quasibinomial))
subc.MIte.results <- pool(subc.MIte.models)
summary_MIte_subc <-summary(subc.MIte.results, conf.int = TRUE)

# IPTW
weighted.MIte.models <- with(data = weighted.MIte.datasets, expr = svyglm(formula1, family = quasibinomial))
weighted.MIte.results <- pool(weighted.MIte.models)
summary_MIte_iptw <-summary(weighted.MIte.results, conf.int = TRUE)

NAFLD.T2D.MIte.HR.est<-c()
NAFLD.T2D.MIte.logHR.est<-c()
NAFLD.T2D.MIte.se<-c()
NAFLD.T2D.MIte.pvalue<-c()
NAFLD.T2D.MIte.lowerCI<-c()
NAFLD.T2D.MIte.upperCI<-c()

NAFLD.T2D.MIte.HR.est[1] <-exp(summary_MIte_greedy[2,2])
NAFLD.T2D.MIte.logHR.est[1]<-summary_MIte_greedy[2,2]
NAFLD.T2D.MIte.se[1]<-summary_MIte_greedy[2,3]
NAFLD.T2D.MIte.pvalue[1]<-summary_MIte_greedy[2,6]
NAFLD.T2D.MIte.lowerCI[1]<- exp(summary_MIte_greedy[2,7])
NAFLD.T2D.MIte.upperCI[1]<- exp(summary_MIte_greedy[2,8])

NAFLD.T2D.MIte.HR.est[2] <-exp(summary_MIte_optimal[2,2])
NAFLD.T2D.MIte.logHR.est[2]<-summary_MIte_optimal[2,2]
NAFLD.T2D.MIte.se[2]<-summary_MIte_optimal[2,3]
NAFLD.T2D.MIte.pvalue[2]<-summary_MIte_optimal[2,6]
NAFLD.T2D.MIte.lowerCI[2]<- exp(summary_MIte_optimal[2,7])
NAFLD.T2D.MIte.upperCI[2]<- exp(summary_MIte_optimal[2,8])

NAFLD.T2D.MIte.HR.est[3] <-exp(summary_MIte_one_two[2,2])
NAFLD.T2D.MIte.logHR.est[3]<-summary_MIte_one_two[2,2]
NAFLD.T2D.MIte.se[3]<-summary_MIte_one_two[2,3]
NAFLD.T2D.MIte.pvalue[3]<-summary_MIte_one_two[2,6]
NAFLD.T2D.MIte.lowerCI[3]<- exp(summary_MIte_one_two[2,7])
NAFLD.T2D.MIte.upperCI[3]<- exp(summary_MIte_one_two[2,8])

NAFLD.T2D.MIte.HR.est[4] <-exp(summary_MIte_full[2,2])
NAFLD.T2D.MIte.logHR.est[4]<-summary_MIte_full[2,2]
NAFLD.T2D.MIte.se[4]<-summary_MIte_full[2,3]
NAFLD.T2D.MIte.pvalue[4]<-summary_MIte_full[2,6]
NAFLD.T2D.MIte.lowerCI[4]<- exp(summary_MIte_full[2,7])
NAFLD.T2D.MIte.upperCI[4]<- exp(summary_MIte_full[2,8])

NAFLD.T2D.MIte.HR.est[5] <-exp(summary_MIte_subc[2,2])
NAFLD.T2D.MIte.logHR.est[5]<-summary_MIte_subc[2,2]
NAFLD.T2D.MIte.se[5]<-summary_MIte_subc[2,3]
NAFLD.T2D.MIte.pvalue[5]<-summary_MIte_subc[2,6]
NAFLD.T2D.MIte.lowerCI[5]<- exp(summary_MIte_subc[2,7])
NAFLD.T2D.MIte.upperCI[5]<- exp(summary_MIte_subc[2,8])

NAFLD.T2D.MIte.HR.est[6] <-exp(summary_MIte_iptw[2,2])
NAFLD.T2D.MIte.logHR.est[6]<-summary_MIte_iptw[2,2]
NAFLD.T2D.MIte.se[6]<-summary_MIte_iptw[2,3]
NAFLD.T2D.MIte.pvalue[6]<-summary_MIte_iptw[2,6]
NAFLD.T2D.MIte.lowerCI[6]<- exp(summary_MIte_iptw[2,7])
NAFLD.T2D.MIte.upperCI[6]<- exp(summary_MIte_iptw[2,8])

NAFLDoutcome.T2D.MIte.result <- 
  structure(list(
    mean  = c(NA, NA, NAFLD.T2D.MIte.HR.est), 
    lower = c(NA, NA, NAFLD.T2D.MIte.lowerCI),
    upper = c(NA, NA, NAFLD.T2D.MIte.upperCI)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

A <- cbind2(c("","Method","Greedy Matching","Optimal Matching","1:2 Matching","Full Matching","Subclassification", "IPTW"),
            c("","Effect", as.character(round(NAFLD.T2D.MIte.HR.est,3))))
B<- cbind2(c("","SE", as.character(round(NAFLD.T2D.MIte.se,3))),
           c("","p-value", ifelse(NAFLD.T2D.MIte.pvalue<0.001, "<0.001", round(NAFLD.T2D.MIte.pvalue,3))))
tabletext <- cbind2(A ,B)

forestplot(tabletext, 
           NAFLDoutcome.T2D.result,new_page = TRUE,
           col = fpColors(box="royalblue",line="darkblue"),
           xlab = "Estimate with 95% CI",
           hrzl_lines= gpar(col = "#444444"),
           is.summary = c(TRUE, TRUE, rep(FALSE, 10)),
           clip = c(0.5,5.5), zero=1,
           txt_gp=fpTxtGp(label = gpar(cex = 0.8),
                          title = gpar(cex = 1),
                          ticks = gpar(cex = 0.7),
                          xlab = gpar(cex = 0.7)),
           grid = structure(seq(1, 5.5, by=0.25), 
                            gp = gpar(lty = 2, col = "#CCCCFF")))
grid.text("Exposure of T2D on NAFLD by MIte", .5, 0.98, gp=gpar(fontface="bold"))

################################
############MIps################
################################

## Greedy matching
greedy.matched.MIps.across.datasets <- matchthem(ps.formula, mi.data,
                                                 approach = 'across',
                                                 method = 'nearest',
                                                 caliper = 0.05)
## Optimal matching
optimal.matched.MIps.across.datasets <- matchthem(ps.formula, mi.data,
                                                  approach = 'across',
                                                  method = 'optimal',
                                                  caliper = 0.05)
## 1:2 matching
one_two.matched.MIps.across.datasets <- matchthem(ps.formula, mi.data,
                                                  approach = 'across',
                                                  method = 'nearest',
                                                  caliper = 0.05, ratio = 2)
## Full matching
full.matched.MIps.across.datasets <- matchthem(ps.formula,mi.data,
                                               approach = 'across',
                                               caliper = 0.05,
                                               method = "full")

## Subclassification with full matching
subc.MIps.across.datasets <- matchthem(ps.formula, mi.data,
                                       approach = 'across',
                                       method = "subclass",
                                       subclass = 5)
subc.MIps.across.datasets
plot(subc.MIps.across.datasets, type = "jitter")

# IPTW
weighted.across.datasets <- weightthem(ps.formula, mi.data,
                                       approach = 'across', 
                                       method = 'ps', estimand = 'ATM')
weighted.across.datasets

bal.tab(greedy.matched.MIps.across.datasets, abs = TRUE)
bal.tab(optimal.matched.MIps.across.datasets, abs = TRUE)
bal.tab(one_two.matched.MIps.across.datasets, abs = TRUE)
bal.tab(full.matched.MIps.across.datasets, abs = TRUE)
bal.tab(subc.MIps.across.datasets, abs = TRUE)
bal.tab(weighted.across.datasets, abs = TRUE)

formula1<-as.formula((paste("NAFLD_status","~ T2D+ ",paste(newvariables,collapse = "+"))))

greedy.matched.MIps.across.models <- with(data = greedy.matched.MIps.across.datasets,
                                          expr = svyglm(formula1, family = quasibinomial))
greedy.matched.MIps.across.results <- pool(greedy.matched.MIps.across.models)
summary_MIps_greedy <- summary(greedy.matched.MIps.across.results, conf.int = TRUE)

optimal.matched.MIps.across.models <- with(data = optimal.matched.MIps.across.datasets,
                                           expr = svyglm(formula1, family =quasibinomial))
optimal.matched.MIps.across.results <- pool(optimal.matched.MIps.across.models)
summary_MIps_optimal <-summary(optimal.matched.MIps.across.results, conf.int = TRUE)

one_two.matched.MIps.across.models <- with(data = one_two.matched.MIps.across.datasets,
                                           expr = svyglm(formula1, family = quasibinomial))
one_two.matched.MIps.across.results <- pool(one_two.matched.MIps.across.models)
summary_MIps_one_two <- summary(one_two.matched.MIps.across.results, conf.int = TRUE)

full.matched.MIps.across.models <- with(data = full.matched.MIps.across.datasets,
                                        expr = svyglm(formula1, family = quasibinomial))
full.matched.MIps.across.results <- pool(full.matched.MIps.across.models)
summary_MIps_full <- summary(full.matched.MIps.across.results, conf.int = TRUE)


subc.MIps.across.models <- with(data = subc.MIps.across.datasets,
                                expr = svyglm(formula1, family = quasibinomial))
subc.MIps.across.results <- pool(subc.MIps.across.models)
summary_MIps_subc <-summary(subc.MIps.across.results, conf.int = TRUE)

# IPTW

weighted.across.models <- with(data = weighted.across.datasets, expr = svyglm(formula1, family = quasibinomial))
weighted.across.results <- pool(weighted.across.models)
summary_MIps_iptw <-summary(weighted.across.results, conf.int = TRUE)

NAFLD.T2D.MIps.HR.est<-c()
NAFLD.T2D.MIps.logHR.est<-c()
NAFLD.T2D.MIps.se<-c()
NAFLD.T2D.MIps.pvalue<-c()
NAFLD.T2D.MIps.lowerCI<-c()
NAFLD.T2D.MIps.upperCI<-c()

NAFLD.T2D.MIps.HR.est[1] <-exp(summary_MIps_greedy[2,2])
NAFLD.T2D.MIps.logHR.est[1]<-summary_MIps_greedy[2,2]
NAFLD.T2D.MIps.se[1]<-summary_MIps_greedy[2,3]
NAFLD.T2D.MIps.pvalue[1]<-summary_MIps_greedy[2,6]
NAFLD.T2D.MIps.lowerCI[1]<- exp(summary_MIps_greedy[2,7])
NAFLD.T2D.MIps.upperCI[1]<- exp(summary_MIps_greedy[2,8])

NAFLD.T2D.MIps.HR.est[2] <-exp(summary_MIps_optimal[2,2])
NAFLD.T2D.MIps.logHR.est[2]<-summary_MIps_optimal[2,2]
NAFLD.T2D.MIps.se[2]<-summary_MIps_optimal[2,3]
NAFLD.T2D.MIps.pvalue[2]<-summary_MIps_optimal[2,6]
NAFLD.T2D.MIps.lowerCI[2]<- exp(summary_MIps_optimal[2,7])
NAFLD.T2D.MIps.upperCI[2]<- exp(summary_MIps_optimal[2,8])

NAFLD.T2D.MIps.HR.est[3] <-exp(summary_MIps_one_two[2,2])
NAFLD.T2D.MIps.logHR.est[3]<-summary_MIps_one_two[2,2]
NAFLD.T2D.MIps.se[3]<-summary_MIps_one_two[2,3]
NAFLD.T2D.MIps.pvalue[3]<-summary_MIps_one_two[2,6]
NAFLD.T2D.MIps.lowerCI[3]<- exp(summary_MIps_one_two[2,7])
NAFLD.T2D.MIps.upperCI[3]<- exp(summary_MIps_one_two[2,8])

NAFLD.T2D.MIps.HR.est[4] <-exp(summary_MIps_full[2,2])
NAFLD.T2D.MIps.logHR.est[4]<-summary_MIps_full[2,2]
NAFLD.T2D.MIps.se[4]<-summary_MIps_full[2,3]
NAFLD.T2D.MIps.pvalue[4]<-summary_MIps_full[2,6]
NAFLD.T2D.MIps.lowerCI[4]<- exp(summary_MIps_full[2,7])
NAFLD.T2D.MIps.upperCI[4]<- exp(summary_MIps_full[2,8])

NAFLD.T2D.MIps.HR.est[5] <-exp(summary_MIps_subc[2,2])
NAFLD.T2D.MIps.logHR.est[5]<-summary_MIps_subc[2,2]
NAFLD.T2D.MIps.se[5]<-summary_MIps_subc[2,3]
NAFLD.T2D.MIps.pvalue[5]<-summary_MIps_subc[2,6]
NAFLD.T2D.MIps.lowerCI[5]<- exp(summary_MIps_subc[2,7])
NAFLD.T2D.MIps.upperCI[5]<- exp(summary_MIps_subc[2,8])

NAFLD.T2D.MIps.HR.est[6] <-exp(summary_MIps_iptw[2,2])
NAFLD.T2D.MIps.logHR.est[6]<-summary_MIps_iptw[2,2]
NAFLD.T2D.MIps.se[6]<-summary_MIps_iptw[2,3]
NAFLD.T2D.MIps.pvalue[6]<-summary_MIps_iptw[2,6]
NAFLD.T2D.MIps.lowerCI[6]<- exp(summary_MIps_iptw[2,7])
NAFLD.T2D.MIps.upperCI[6]<- exp(summary_MIps_iptw[2,8])

NAFLDoutcome.T2D.MIps.result <- 
  structure(list(
    mean  = c(NA, NA, NAFLD.T2D.MIps.HR.est), 
    lower = c(NA, NA, NAFLD.T2D.MIps.lowerCI),
    upper = c(NA, NA, NAFLD.T2D.MIps.upperCI)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

tabletext <- cbind(c("","Method","Greedy Matching","Optimal Matching","1:2 Matching","Full Matching","Subclassification", "IPTW"), c("","Effect", as.character(round(NAFLD.T2D.MIps.HR.est,3))), c("","SE", as.character(round(NAFLD.T2D.MIps.se,3))),c("","p-value", ifelse(NAFLD.T2D.MIps.pvalue<0.001, "<0.001", round(NAFLD.T2D.MIps.pvalue,3))))

forestplot(tabletext, 
           NAFLDoutcome.T2D.result,new_page = TRUE,
           col = fpColors(box="royalblue",line="darkblue"),
           xlab = "Estimate with 95% CI",
           hrzl_lines= gpar(col = "#444444"),
           is.summary = c(TRUE, TRUE, rep(FALSE, 10)),
           clip = c(0.5,5.5), zero=1,
           txt_gp=fpTxtGp(label = gpar(cex = 0.8),
                          title = gpar(cex = 1),
                          ticks = gpar(cex = 0.8),
                          xlab = gpar(cex = 0.7)),
           grid = structure(seq(1, 5.5, by=0.25), 
                            gp = gpar(lty = 2, col = "#CCCCFF")))
grid.text("Exposure of T2D on NAFLD by MIps", .5, 0.98, gp=gpar(fontface="bold"))

NAFLD.T2D<-data.frame(method=c(rep("CCA",6),rep("MIte",6),rep("MIps",6)),PSmethod=c(rep(c("Greedy Matching","Optimal Matching","1:2 Matching","Full Matching","Subclassification", "IPTW"),3)),est=c(round(NAFLD.T2D.HR.est,3),round(NAFLD.T2D.MIte.HR.est,3),round(NAFLD.T2D.MIps.HR.est,3)),log.est=c(round(NAFLD.T2D.logHR.est,3),round(NAFLD.T2D.MIte.logHR.est,3),round(NAFLD.T2D.MIps.logHR.est,3)),se=c(round(NAFLD.T2D.se,3),round(NAFLD.T2D.MIte.se,3),round(NAFLD.T2D.MIps.se,3)),pvalue=c(round(NAFLD.T2D.pvalue,3),round(NAFLD.T2D.MIte.pvalue,3),round(NAFLD.T2D.MIps.pvalue,3)),lowerCI=c(round(NAFLD.T2D.lowerCI,3),round(NAFLD.T2D.MIte.lowerCI,3),round(NAFLD.T2D.MIps.lowerCI,3)),upperCI=c(round(NAFLD.T2D.upperCI,3),round(NAFLD.T2D.MIte.upperCI,3),round(NAFLD.T2D.MIps.upperCI,3)))

NAFLD.T2D.greedy<-subset(NAFLD.T2D,PSmethod=="Greedy Matching")
NAFLD.T2D.optimal<-subset(NAFLD.T2D,PSmethod=="Optimal Matching")
NAFLD.T2D.one_two<-subset(NAFLD.T2D,PSmethod=="1:2 Matching")
NAFLD.T2D.full<-subset(NAFLD.T2D,PSmethod=="Full Matching")
NAFLD.T2D.subc<-subset(NAFLD.T2D,PSmethod=="Subclassification")
NAFLD.T2D.iptw<-subset(NAFLD.T2D,PSmethod=="IPTW")


Forstplot<-function(method,title){
  A <- cbind2(c("","Method","CCA","MIte","MIps"), 
                     c("","Effect", as.character(method$est))) 
  B <-cbind2(c("","SE", as.character(method$se)),
                     c("","p-value", ifelse(method$pvalue<0.001, "<0.001", method$pvalue)))
tabletext <- cbind2(A,B)
  NAFLD.T2D.result <- 
    structure(list(
      mean  = c(NA, NA, method$est), 
      lower = c(NA, NA, method$lowerCI),
      upper = c(NA, NA, method$upperCI)),
      .Names = c("mean", "lower", "upper"), 
      row.names = c(NA, -5L), 
      class = "data.frame")
  
  forestplot(tabletext, 
             NAFLD.T2D.result,
             new_page = TRUE,
             xlab = "Estimate with 95% CI",
             col = fpColors(box="royalblue",line="darkblue"),
             hrzl_lines= gpar(col = "#444444"),
             is.summary = c(TRUE, TRUE, rep(FALSE, 18)),
             clip = c(0.5,5.5), zero=1,
             txt_gp=fpTxtGp(label = gpar(cex = 0.8),
                            title = gpar(cex = 1),
                            ticks = gpar(cex = 0.7),
                            xlab = gpar(cex = 0.7)),
             grid = structure(seq(1, 5.5, by=0.5), 
                              gp = gpar(lty = 2, col = "#CCCCFF")),
             title=title)
}

Forstplot(NAFLD.T2D.greedy,"Exposure of T2D on of Greedy Matching")
Forstplot(NAFLD.T2D.optimal,"Exposure of T2D on of Optimal Matching")
Forstplot(NAFLD.T2D.one_two,"Exposure of T2D on of 1:2 Matching")
Forstplot(NAFLD.T2D.full,"Exposure of T2D on of Full Matching")
Forstplot(NAFLD.T2D.subc,"Exposure of T2D on of Subclassification")
Forstplot(NAFLD.T2D.iptw,"Exposure of T2D on of IPTW ")
