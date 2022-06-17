library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

load(file = here::here("output","cohort_long.RData"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

## Count variables for demographic tables
dt.tv[,postVTE90.perepisode := max(post.VTE & end <=90,na.rm = T), keyby = .(patient_id, end.fu)]
dt.tv[!is.finite(postVTE90.perepisode),postVTE90.perepisode := 0]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[ ,postCOVID30.perepisode := event == 1 & end <=30]

dt.tv[,postCOVID30.perepisode := max(postCOVID30.perepisode,na.rm = T), keyby = .(patient_id, end.fu)]
dt.tv[!is.finite(postCOVID30.perepisode),postCOVID30.perepisode := 0]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv[,event.los := (died == T)*2]
dt.tv[event.los !=2,event.los := discharged == T]
dt.tv[, postop.los.cohort := start>=0 & tstop <= los.end & end <= 90 & any.op ==T]
dt.tv[,age.cat := factor(age.cat, order = F)]

n.type.events <- sort(unique(dt.tv[(postop.los.cohort) ,event.los]))[-1]

post.op.LOS.model <-  flexsurv::flexsurvreg(survival::Surv(start,end, event.los == 1) ~ Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + postcovid + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, 
                                           data = dt.tv[(postop.los.cohort)],
                                           dist = 'gengamma')

#  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end, event.los == i) ~ Abdominal + Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + postcovid + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,data = dt.tv[(postop.los.cohort)], model = T))


new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(30,8*length(procedures)),
                                                'event' = rep(F,8*length(procedures)),
                                                'Abdominal' = c(rep(T,8),rep(F,40)),
                                                'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
                                                'Obstetrics'=c(rep(F,16),rep(T,8),rep(F,24)),
                                                'Orthopaedic'=c(rep(F,24),rep(T,8),rep(F,16)),
                                                'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
                                                'Vascular'=c(rep(F,40),rep(T,8)),
                                                'postcovid' = as.numeric(rep(c(rep(F,4),rep(T,4)), times = length(procedures))),
                                                'age.cat' = rep('(50,70]',8*length(procedures)),
                                                'sex' = rep('F',8*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                'region' = rep("East Midlands",8*length(procedures)),
                                                'Current.Cancer' = rep(T,8*length(procedures)),
                                                'Emergency' =   rep(F,8*length(procedures)),
                                                'Charl12' =  rep('Single',8*length(procedures)),
                                                'recentCOVID' = rep(F,8*length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))



mean.adjusted.los <- predict(post.op.LOS.model, newdata = new.data.postop.covid, type = "response", conf.int = T)
 # matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv', model = 'post.op.LOS.model', newdata = 'new.data.postop.covid', day = 14), byrow = T, ncol = 4)

mean.adjusted.los <- matrix(apply(mean.adjusted.los,1,function(x) paste0(round(x[1], digits = 1), 
                                                                  " days, (",
                                                                  round(x[2], digits = 1),
                                                                  ',',
                                                                  round(x[3], digits = 1),')')), byrow = T, ncol = 4)

colnames(mean.adjusted.los) <- paste0('Wave_',1:4)
rownames(mean.adjusted.los) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))
save(mean.adjusted.los,post.op.LOS.model, file = here::here("output","postoplos.RData"))
