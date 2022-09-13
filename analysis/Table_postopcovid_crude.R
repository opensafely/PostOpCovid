library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'sex','age.cat','imd5','wave','bmi.cat',
                'vaccination.status.factor','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

summary(dt.tv)

#Could be removed after data.manage rerun
dt.tv[, region:= as.factor(region)]
#bin.cov <- c(procedures,'sex','Current.Cancer','Emergency','recentCOVID','previousCOVID')
#dt.tv[,(bin.cov) := lapply(.SD, function(x) data.table::fifelse(is.na(x),F,x)), .SDcols = c(bin.cov)]

crude.covid.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[(postop.covid.cohort),
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[(postop.covid.cohort),as.factor(get(covariates[i]))]),
                                                                                        cuminc.km(covariates[i], niter = 2)[,2:5])))

names(crude.covid.cov) <- c("Characteristic","Level","Number at risk",
                            "Number of events","30 day Cumulative Risk adjusted for censoring","30 day Cumulative Risk adjusted for death and emergency readmission")

data.table::fwrite(crude.covid.cov, file = here::here("output","postopcovid_crude.csv"))
save(crude.covid.cov, file = here::here("output","postopcovid_crude.RData"))
