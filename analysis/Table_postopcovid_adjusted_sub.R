library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(feather::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')


covariates <- c(procedures,'sex','age.cat','bmi.cat','imd5','wave',
                'vaccination.status.factor','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID','region')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
      (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]

post.op.covid.model.waves.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Colectomy + Cholecystectomy + HipReplacement + KneeReplacement + age.cat + sex + bmi.cat + imd5 + vaccination.status.factor + region + Current.Cancer + Emergency*wave + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort) & sub.op == T], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.waves.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelwavessub.csv"))


new.data.postop <- data.table::data.table(
  'start' = rep(0,8*length(procedures)),
  'end' = rep(30,8*length(procedures)),
  'event' = rep(F,8*length(procedures)),
  'Colectomy' = c(rep(T,8),rep(F,24)),
  'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
  'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
  'KneeReplacement'=c(rep(F,24),rep(T,8)),
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

cuminc.adjusted.waves.sub <- 
  matrix(cuminc.cox(n.type.events = n.type.events,
                    dt = 'dt.tv[(postop.covid.cohort) & sub.op == T]',
                    model = 'post.op.covid.model.waves.sub', 
                    newdata = 'new.data.postop',
                    day = 30), byrow = T, ncol = 4)

colnames(cuminc.adjusted.waves.sub) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.waves.sub) <- paste0(c('Elective','Emergency'),rep(procedures, each = 2))

data.table::fwrite(cuminc.adjusted.waves.sub, file = here::here("output","postopcovid_adjusted_waves_sub.csv"))

#############################################################################################

data.table::setkey(dt.tv,"patient_id","tstart","tstop")
post.op.covid.model.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Colectomy + Cholecystectomy + HipReplacement + KneeReplacement +
                                                      age.cat + sex + bmi.cat + imd5 + 
                                                      vaccination.status.factor + region + Current.Cancer + 
                                                      Emergency + wave + Charl12 + recentCOVID + previousCOVID, 
                                                    id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort)  & sub.op == T], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelsub.csv"))


adjusted.cuminc.sub <-  data.table::as.data.table(foreach::foreach(predi = 1:length(covariates), .combine = 'c', .inorder = T) %do% {
                           newdata.rows <- length(unique(dt.tv[!is.na(get(covariates[predi])) ,get(covariates[predi])]))
                           
   
                           newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                                                                  'end' = rep(30,newdata.rows),
                                                                  'event' = rep(F,newdata.rows),
                                                                  'patient_id' = 1:newdata.rows,
                                                                  'Colectomy' = c(rep(T,newdata.rows)),
                                                                  'Cholecystectomy'=c(rep(F,newdata.rows)),
                                                                  'HipReplacement'=c(rep(F,newdata.rows)),
                                                                  'KneeReplacement'=c(rep(F,newdata.rows)),
                                                                  'age.cat' = rep('(50,70]',newdata.rows),
                                                                  'sex' = rep('F',newdata.rows),
                                                                  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],newdata.rows),
                                                                  'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows),
                                                                  'wave' = rep(paste0('Wave_',4),times = newdata.rows),
                                                                  'vaccination.status.factor' = rep('3',newdata.rows),
                                                                  'region' = rep("East Midlands",newdata.rows),
                                                                  'Current.Cancer' = rep(T,newdata.rows),
                                                                  'Emergency' =  rep(F,newdata.rows),
                                                                  'Charl12' =  rep('Single',newdata.rows),
                                                                  'recentCOVID' = rep(F,newdata.rows),
                                                                  'previousCOVID' = rep(F,newdata.rows)
                                                                  )
                            if ( predi <= length(procedures)) {
                                newdata.pred[,(procedures) := F]
                                newdata.pred[,(procedures[predi]) := c(F,T)]
                            } else {

                            
                          #  newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                          #                                         'end' = rep(30,newdata.rows),
                          #                                         'event' = rep(F,newdata.rows),
                          #                                         'patient_id' = 1:newdata.rows)
                          #  if ( predi > length(procedures)) {
                          #    newdata.pred[,(procedures) := lapply(procedures, function(x) x == procedures[which.max(dt.tv[,lapply(.SD,sum,na.rm = T), .SDcols = c(procedures)])])] } else {
                          #      newdata.pred[,(procedures) := lapply(procedures, function(x) x == covariates[predi] & patient_id > 1)]
                          #    }
                           
                          #  newdata.pred[,(covariates[-c(1:length(procedures))]) := lapply(((length(procedures)+1):length(covariates)), function(i.c) {
                          #    if(is.factor(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                          #      as.character(rep(max.category(i.c),newdata.rows))
                          #    } else if(is.logical(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                          #      as.logical(rep(max.category(i.c),newdata.rows))
                          #    } else if(is.numeric(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
                          #      is.numeric(rep(max.category(i.c),newdata.rows))
                          #    } else {
                          #      rep(max.category(i.c),newdata.rows)
                          #    }
                          #  })]
                           
                          #  #names(newdata.pred) <- c('start','end','event', covariates,'patient_id') 
                           if(is.factor(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  as.character(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.logical(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  as.logical(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.numeric(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  is.numeric(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else {
                              newdata.pred[,(covariates[predi]) := sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T))]
                            }
                            }
                           cuminc.cox(n.type.events = n.type.events,
                                      dt = 'dt.tv[(postop.covid.cohort) & sub.op == T]', 
                                      model = 'post.op.covid.model.sub', 
                                      newdata = 'newdata.pred',
                                      day = 30)
                         })


save(cuminc.adjusted.waves.sub,post.op.covid.model.sub,adjusted.cuminc.sub, file = here::here("output","postopcovid_adjusted_sub.RData"))
data.table::fwrite(cuminc.adjusted.sub, file = here::here("output","postopcovid_adjusted_sub.csv"))
