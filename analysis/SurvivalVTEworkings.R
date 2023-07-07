library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

##########

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave','Major.op','vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

summary(dt.tv)

dt.tv[, region:= as.factor(region)]

###Calculate crude VTE rates
dt.tv[, vte := max(as.numeric(event.VTE == 1),na.rm = TRUE),keyby= patient_id]

# create column 'vte' indicating if a patient had VTE

#dt.tv[, .N, by = list( Major.op, vte)][,.(Total = sum (N),Events = sum(N[vte==1])), by = Major.op][,.( covariate= 'Major.op', levels= Major.op,Total,Events, Crude.Rate = 100 * Events / Total)]

crude.vte.rates<-cbind(rbindlist(lapply(covariates,function(x)dt.tv[,head(.SD,1),keyby=.(patient_id, end.fu)][,
                                   
                                .N,
                                keyby = c( x, 'vte')][,
                                                           .(Total = sum (N),
                                                             Events = sum(N[vte==1])), 
                                                            keyby = x][,.( covariate= x, 
                                                                               levels= get(x),
                                                                               Total,
                                                                               Events, 
                                                                               Crude.Rate = 100 * Events / Total)] )),



rbindlist(lapply(covariates,function(x)cbind((data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event.VTE==1) ~ get(x), 
                                                data = dt.tv[(postcovid.VTE.cohort) & !is.na(get(x))  & start >= 0,],
                                                id = patient_id), 
                              times = 0,
                              extend = T)[c('n.risk')])),
(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event.VTE==1) ~ get(x), 
                                                data = dt.tv[(postcovid.VTE.cohort) & !is.na(get(x)) & start >= 0,],
                                                id = patient_id), 
                              times = 0,
                              extend = T)[c('n.event')])),
100*round(1 - (data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event.VTE==1) ~ get(x), 
                                                           data = dt.tv[(postcovid.VTE.cohort) & !is.na(get(x)) ,],
                                                           id = patient_id), times = 30,
                                         extend = T)[c('surv')])))))))



data.table::fwrite(crude.vte.rates, file = here::here("output","crude_vte_rates.csv"))
save(crude.vte.rates, file = here::here("output","crude_vte_rates.RData"))