library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

load(file = here::here("output","cohort_long.RData"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,age.cat := factor(age.cat, order = F)]

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

crude.covid.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[(postop.covid.cohort),
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[(postop.covid.cohort),as.factor(get(covariates[i]))]),
                                                                                        cuminc.km(covariates[i], niter = 2)[,2:5])))

names(crude.covid.cov) <- c("Characteristic","Level","Number at risk",
                            "Number of events","30 day Cumulative Risk adjusted for censoring","30 day Cumulative Risk adjusted for death and emergency readmission")

save(crude.covid.cov, file = here::here("output","postopcovid_crude.RData"))