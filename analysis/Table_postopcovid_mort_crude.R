library(foreach)
library(data.table)
library(magrittr)
library(survival)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','postcovid','wave','Major.op',
                'vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

summary(dt.tv)

#Could be removed after data.manage rerun
dt.tv[, region:= as.factor(region)]
dt.tv[,postcovid := postcovid == 1]
#bin.cov <- c(procedures,'sex','Current.Cancer','Emergency','recentCOVID','previousCOVID')
#dt.tv[,(bin.cov) := lapply(.SD, function(x) data.table::fifelse(is.na(x),F,x)), .SDcols = c(bin.cov)]

crude.HR <- function(x) {vapply(1:length(coef(x)), FUN.VALUE = '', function(i) paste0(round(coef(x)[i],2),' (',paste(round(confint(x)[i,],2), collapse = ','),')'))}
crude.mort.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[start >=0 & any.op == T,
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[start >=0 & any.op == T,as.factor(get(covariates[i]))]),
                                                                                        cuminc.km.mort(covariates[i], niter = 2)[,2:4],
                                                                                       c('Ref',coxph(formula = as.formula(paste0('Surv(start,end,died) ~ ',covariates[i])),  
                                                                                                     id = patient_id,
                                                                                                     data = dt.tv[start >=0 & any.op == T]) %>% crude.HR()))))

names(crude.mort.cov) <- c("Characteristic",
                           "Level",
                           "Number at risk",
                            "Number of events",
                           "30 day Cumulative Risk adjusted for censoring",
                           "Crude Hazard ratio (95% CI)")

data.table::fwrite(crude.mort.cov, file = here::here("output","postopmort_crude.csv"))
save(crude.mort.cov, file = here::here("output","postopmort_crude.RData"))
