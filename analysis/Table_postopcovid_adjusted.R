library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

load(file = here::here("output","cohort_long.RData"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

crude.covid.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[(postop.covid.cohort),
                                                                                                                              as.factor(get(covariates[i]))]))), 
                                                                                        levels(dt.tv[(postop.covid.cohort),as.factor(get(covariates[i]))]),
                                                                                        cuminc.km(covariates[i], niter = 2)[,2:5])))

names(crude.covid.cov) <- c("Characteristic","Level","Number at risk",
                            "Number of events","30 day Cumulative Risk adjusted for censoring","30 day Cumulative Risk adjusted for death and emergency readmission")

data.table::setkey(dt.tv,"patient_id","tstart","tstop")

dt.tv[,age.cat := factor(age.cat, order = F)]
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

post.op.covid.model <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Abdominal + Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort)], model = T))

new.data.postop <- data.table::data.table(
  'start' = rep(0,8*length(procedures)),
  'end' = rep(30,8*length(procedures)),
  'event' = rep(F,8*length(procedures)),
  'Abdominal' = c(rep(T,8),rep(F,40)),
  'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
  'Obstetrics'=c(rep(F,16),rep(T,8),rep(F,24)),
  'Orthopaedic'=c(rep(F,24),rep(T,8),rep(F,16)),
  'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
  'Vascular'=c(rep(F,40),rep(T,8)),
  'age.cat' = rep('(50,70]',8*length(procedures)),
  'sex' = rep('F',8*length(procedures)),
  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
  'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
  'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
  'vaccination.status.factor' = rep('3',8*length(procedures)),
  'region' = rep("East Midlands",8*length(procedures)),
  'Current.Cancer' = rep(T,8*length(procedures)),
  'Emergency' =  rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
  'Charl12' =  rep('Single',8*length(procedures)),
  'recentCOVID' = rep(F,8*length(procedures)),
  'previousCOVID' = rep(F,8*length(procedures)),
  'patient_id' = 1:(8*length(procedures)))


cuminc.adjusted.waves <- 
  matrix(cuminc.cox(n.type.events = n.type.events,
                    dt = 'dt.tv',
                    model = 'post.op.covid.model', 
                    newdata = 'new.data.postop',
                    day = 30), byrow = T, ncol = 4)

colnames(cuminc.adjusted.waves) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.waves) <- paste0(c('Elective','Emergency'),rep(procedures, each = 2))


data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.covid.model <- 
  lapply(1:3, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                          data = dt.tv[(postop.covid.cohort)], model = T))


n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

adjusted.cuminc <-  foreach::foreach(predi = 1:length(covariates), .combine = 'c', .inorder = T) %do% {
                           newdata.rows <- length(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])]))
                           
                           newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                                                                  'end' = rep(30,newdata.rows),
                                                                  'event' = rep(F,newdata.rows),
                                                                  'patient_id' = 1:newdata.rows)
                           if ( predi > length(procedures)) {
                             newdata.pred[,(procedures) := lapply(procedures, function(x) x == procedures[which.max(dt.tv[,lapply(.SD,sum,na.rm = T), .SDcols = c(procedures)])])] } else {
                               newdata.pred[,(procedures) := lapply(procedures, function(x) x == covariates[predi] & patient_id > 1)]
                             }
                           
                           newdata.pred[,(covariates[-c(1:length(procedures))]) := lapply(((length(procedures)+1):length(covariates)), function(i.c) {
                             if(is.factor(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                               as.character(rep(max.category(i.c),newdata.rows))
                             } else if(is.logical(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                               as.logical(rep(max.category(i.c),newdata.rows))
                             } else if(is.numeric(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                               is.numeric(rep(max.category(i.c),newdata.rows))
                             } else {
                               rep(max.category(i.c),newdata.rows)
                             }
                           })]
                           
                           #names(newdata.pred) <- c('start','end','event', covariates,'patient_id') 
                           if(is.factor(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                             newdata.pred[,(covariates[predi]) :=  as.character(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                           } else if(is.logical(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                             newdata.pred[,(covariates[predi]) :=  as.logical(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                           } else if(is.numeric(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                             newdata.pred[,(covariates[predi]) :=  is.numeric(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                           } else {
                             newdata.pred[,(covariates[predi]) := sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T))]
                           }
                           
                           cuminc.cox(n.type.events = n.type.events,
                                      dt = 'dt.tv', 
                                      model = 'post.op.covid.model', 
                                      newdata = 'newdata.pred',
                                      day = 30)
                         }


save(cuminc.adjusted.waves,post.op.covid.model,adjusted.cuminc, file = here::here("output","postopcovid_adjusted.RData"))