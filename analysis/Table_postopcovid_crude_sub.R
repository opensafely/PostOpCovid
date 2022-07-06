library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')


covariates <- c(procedures,'sex','age.cat','bmi.cat','imd5','wave',
                'vaccination.status.factor','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID','region')

dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]


data.table::setkey(dt.tv,patient_id,tstart,tstop)

crude.covid.cov.sub <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[(postop.covid.cohort),
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[(postop.covid.cohort),as.factor(get(covariates[i]))]),
                                                                                        cuminc.km.sub(covariates[i], niter = 2)[,2:5])))

names(crude.covid.cov.sub) <- c("Characteristic","Level","Number at risk",
                            "Number of events","30 day Cumulative Risk adjusted for censoring","30 day Cumulative Risk adjusted for death and emergency readmission")

data.table::fwrite(crude.covid.cov.sub, file = here::here("output","postopcovid_crude_sub.csv"))

save(crude.covid.cov.sub, file = here::here("output","postopcovid_crude_sub.RData"))
