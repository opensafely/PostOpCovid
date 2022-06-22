library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

load(file = here::here("output","cohort_long.RData"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv',
           ID = 'patient_id',
           start.DTTM = 'tstart',
           group = 'c("patient_id","end.fu")',
           var.cols = paste0('c("discharge.date"',paste(covariates,collapse = '","'),')')


dt.tv[, `:=`(start = tstart - discharge.date,
                    end = tstop - discharge.date)]

dt.tv <- dt.tv[start>=0 & end <=90]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,final.date := readmit.end]
dt.tv[is.finite(end.fu) & end.fu < final.date, final.date := end.fu]
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date',aggregate.cols = 'final.date',id.vars = c("patient_id","end.fu"))


dt.tv[,event.readmit :=0]
dt.tv[emergency_readmitdate  == tstop & COVIDpositivedate != tstop, event.readmit := 1]
dt.tv[emergency_readmitdate  == tstop & COVIDpositivedate == tstop, event.readmit := 2]
dt.tv[date_death_ons == tstop & event.readmit != 1, event.readmit := 3]

dt.tv[, postop.readmit.cohort := start>=0 & tstop <= final.date & end <= 90]

dt.tv[(postop.readmit.cohort) & start ==0 ,any.op.readmit := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op.readmit), any.op.readmit := F]
dt.tv[, any.op.readmit := any.op.readmit > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op.readmit := cummax(any.op.readmit), keyby = .(patient_id, end.fu)]

dt.tv[, postop.readmit.cohort := start>=0 & tstop <= final.date & end <= 90 & any.op.readmit == T]

data.table::setkey(dt.tv, patient_id, end.fu, start)

post.op.readmit.model <- 
  lapply(1:3, function(i) survival::coxph(survival::Surv(start,end,event.readmit==i) ~ Abdominal + Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + postcovid*wave + age.cat + sex + bmi.cat + imd5 + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                          data = dt.tv[(postcovid.readmit.cohort)], model = T))

data.table::fwrite(broom::tidy(post.op.readmit.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopreadmitmodel.csv"))


n.type.events <- sort(unique(dt.tv[(postop.readmit.cohort) ,event]))[-1]

new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(30,8*length(procedures)),
                                                'event' = rep(F,8*length(procedures)),
                                                'Abdominal' = c(rep(T,8),rep(F,40)),
                                                'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
                                                'Obstetrics'=c(rep(F,16),rep(T,8),rep(F,24)),
                                                'Orthopaedic'=c(rep(F,24),rep(T,8),rep(F,16)),
                                                'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
                                                'Vascular'=c(rep(F,40),rep(T,8)),
                                                'postcovid' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
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

cuminc.adjusted.readmit <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.readmit.cohort)]', model = 'post.op.readmit.model', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.readmit) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.readmit) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.readmit.model,cuminc.adjusted.readmit, file = here::here("output","postopreadmit.RData"))