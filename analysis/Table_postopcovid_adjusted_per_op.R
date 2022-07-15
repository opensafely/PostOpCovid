library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

procedures <- c('Cardiac')

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

data.table::setkey(dt.tv,patient_id,tstart,tstop)


covariates <- c('sex','age.cat','bmi.cat','imd5','wave',
                'vaccination.status.factor','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID','region')

data.table::setkey(dt.tv,patient_id,tstart,tstop)


n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

post.op.covid.model.per.op <- 
  lapply(1:length(procedures), function(p) lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ 
                                                                                               age.cat + sex + bmi.cat + imd5 + 
                                                                                               vaccination.status.factor + region + 
                                                                                               Current.Cancer + Emergency + wave + 
                                                                                               Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort) & (get(procedures[p])==T)], model = T)))


#############################################################################################

data.table::setkey(dt.tv,"patient_id","tstart","tstop")


adjusted.cuminc.per.op <-  data.table::as.data.table(foreach::foreach(proc = 1:length(procedures)) %:% foreach::foreach(predi = 1:length(covariates), .combine = 'rbind', .inorder = T) %do% {
                           newdata.rows <- length(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])]))
                           
   
                           newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                                                                  'end' = rep(30,newdata.rows),
                                                                  'event' = rep(F,newdata.rows),
                                                                  'patient_id' = 1:newdata.rows,
                                                                  'age.cat' = rep(levels(dt.tv$age.cat)[1],newdata.rows),
                                                                  'sex' = rep('F',newdata.rows),
                                                                  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[1],newdata.rows),
                                                                  'imd5' = rep(levels(dt.tv$imd5)[1], newdata.rows),
                                                                  'wave' = rep(paste0('Wave_',1),times = newdata.rows),
                                                                  'vaccination.status.factor' = rep('0',newdata.rows),
                                                                  'region' = rep(levels(dt.tv$region)[1],newdata.rows),
                                                                  'Current.Cancer' = rep(F,newdata.rows),
                                                                  'Emergency' =  rep(F,newdata.rows),
                                                                  'Charl12' =  rep(levels(dt.tv$Charl12)[1],newdata.rows),
                                                                  'recentCOVID' = rep(F,newdata.rows),
                                                                  'previousCOVID' = rep(F,newdata.rows)
                                                                  )


                            
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
                            
                           
                           
                           death.risk.30day <- predict(object = post.op.covid.model.per.op[[proc]][[3]], 
                                                       newdata = newdata.pred, type = 'expected',se.fit = T)
                           
                           readmit.risk.30day <- predict(object = post.op.covid.model.per.op[[proc]][[2]], 
                                                       newdata = newdata.pred, type = 'expected',se.fit = T)
                           
                           covid.risk.30day <- predict(object = post.op.covid.model.per.op[[proc]][[1]], 
                                   newdata = newdata.pred, type = 'expected',se.fit = T)
                           
                           cbind(matrix(paste0(round((1- exp(-covid.risk.30day$fit))*100,3),
                                               ' (', round((1 - exp(-(covid.risk.30day$fit - 1.96*covid.risk.30day$se.fit)))*100,3),',',
                                               round((1 - exp(-(covid.risk.30day$fit + 1.96*covid.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows),
                           cuminc.cox(n.type.events = n.type.events,
                                      dt = paste0('dt.tv[(postop.covid.cohort) & (get(procedures[[',proc,']])==T)]'), 
                                      model = paste0('post.op.covid.model.per.op[[',proc,']]'), 
                                      newdata = 'newdata.pred',
                                      day = 30),
                           matrix(paste0(round((1- exp(-readmit.risk.30day$fit))*100,3),
                                         ' (', round((1 - exp(-(readmit.risk.30day$fit - 1.96*readmit.risk.30day$se.fit)))*100,3),',',
                                         round((1 - exp(-(readmit.risk.30day$fit + 1.96*readmit.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows),
                           matrix(paste0(round((1- exp(-death.risk.30day$fit))*100,3),
                                         ' (', round((1 - exp(-(death.risk.30day$fit - 1.96*death.risk.30day$se.fit)))*100,3),',',
                                         round((1 - exp(-(death.risk.30day$fit + 1.96*death.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows)
                           )
                         })


save(post.op.covid.model.per.op,adjusted.cuminc.per.op, file = here::here("output","postopcovid_adjusted_per_op.RData"))

#lapply(1:length(procedures), function(p) data.table::fwrite(adjusted.cuminc.per.op[[p]], file = here::here("output",paste0("postopcovid_adjusted_per_op",procedures[p],".csv"))))
data.table::fwrite(adjusted.cuminc.per.op, file = here::here("output","postopcovid_adjusted_per_op.csv"))
