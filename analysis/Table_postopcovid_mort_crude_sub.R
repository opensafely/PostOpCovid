library(foreach)
library(data.table)
library(magrittr)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures.sub <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')

dt.tv[,(procedures.sub) := lapply(.SD,function(x) x==1), .SDcols = (procedures.sub)]

covariates <- c(procedures.sub,'age.cat','sex','bmi.cat','imd5','postcovid','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

dt.tv[,postcovid := postcovid == 1]

dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)  ]


data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))
data.table::setkey(dt.tv,patient_id,tstart,tstop)
library(survival)
crude.HR <- function(x) {vapply(1:length(coef(x)), FUN.VALUE = '', function(i) paste0(round(coef(x)[i],2),' (',paste(round(confint(x)[i,],2), collapse = ','),')'))}
crude.mort.cov.sub <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[start >=0 & sub.op == T ,
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[start >=0 & sub.op == T ,as.factor(get(covariates[i]))]),
                                                                                        cuminc.km.sub(covariates[i], niter = 2)[,2:4],
                                                                                        c('Ref',coxph(formula = as.formula(paste0('Surv(start,end,died) ~ ',covariates[i])),  
                                                                                              id = patient_id,
                                                                                              data = dt.tv[start >=0 & sub.op == T ]) %>% crude.HR()))))

names(crude.mort.cov.sub) <- c("Characteristic",
                               "Level",
                               "Number at risk",
                            "Number of events",
                            "30 day Cumulative Risk adjusted for censoring",
                            "Crude Hazard ratio (95% CI)")

data.table::fwrite(crude.mort.cov.sub, file = here::here("output","postopmort_crude_sub.csv"))

save(crude.mort.cov.sub, file = here::here("output","postopmort_crude_sub.RData"))
