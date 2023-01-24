library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')


covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','postcovid','wave','Major.op',
                'vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','postop.covid.cohort','end.fu','any.op','died'))]

dt.tv[,(drop.vars) := NULL]
dt.tv <- dt.tv[start >=0 & any.op == T]

dt.tv[!is.finite(Major.op), Major.op := F]
dt.tv[,postcovid := postcovid == 1]

gc()
n.type.events <- 1 #sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]


#############################################################################################

data.table::setkey(dt.tv,"patient_id","tstart","tstop")
post.op.mort.model <- list()
post.op.mort.model[[1]] <- survival::coxph(survival::Surv(start,end,died) ~ Abdominal  + 
                                                      Obstetrics + CardioThoracicVascular + 
                                                      age.cat + sex  + bmi.cat + imd5 + postcovid +  wave +  Major.op + 
                                                      vaccination.status.factor  + region +  Current.Cancer + 
                                                      Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                                    id = patient_id,
                                                    data = dt.tv[start >=0 & any.op == T], model = T)


adjusted.mort.cuminc <-  data.table::as.data.table(foreach::foreach(predi = 1:length(covariates),
                                                               .combine = 'rbind', 
                                                               .inorder = T) %do% {
                           newdata.rows <- length(unique(dt.tv[start >=0 & any.op == T & !is.na(get(covariates[predi])),get(covariates[predi])]))
                           
   
                           newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                                                                  'end' = rep(30,newdata.rows),
                                                                  'died' = rep(F,newdata.rows),
                                                                  'patient_id' = 1:newdata.rows,
                                                                  'Abdominal' = rep(T,newdata.rows),
                                                               #   'Cardiac'= rep(F,newdata.rows),
                                                                  'Obstetrics'=rep(F,newdata.rows),
                                                                  'Orthopaedic'=rep(F,newdata.rows),
                                                                #  'Thoracic'=rep(F,newdata.rows),
                                                                  'CardioThoracicVascular'=rep(F,newdata.rows),
                                                                  'age.cat' = rep('(50,70]',newdata.rows),
                                                                  'sex' = rep('F',newdata.rows),
                                                                  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],newdata.rows),
                                                                  'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows),
                                                                  'postcovid' = rep(F,newdata.rows),
                                                                  'wave' = rep(paste0('Wave_',3),times = newdata.rows),
                                                                  'Major.op' = rep(T,newdata.rows),
                                                                  'vaccination.status.factor' = rep('3',newdata.rows),
                                                                  'region' = rep("East Midlands",newdata.rows),
                                                                  'Current.Cancer' = rep(T,newdata.rows),
                                                                  'LOS.bin' = rep(F,newdata.rows),
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
                              newdata.pred[,(covariates[predi]) :=  as.character(sort(unique(dt.tv[ start >=0 & any.op == T & !is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.logical(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  as.logical(sort(unique(dt.tv[start >=0 & any.op == T & !is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.numeric(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  is.numeric(sort(unique(dt.tv[start >=0 & any.op == T & !is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else {
                              newdata.pred[,(covariates[predi]) := sort(unique(dt.tv[start >=0 & any.op == T & !is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T))]
                            }
                            }
                           
                           
                          #  samples <- foreach::foreach(i = 1:1000, .combine = cbind, .multicombine = T, .inorder = F, .verbose = F,
                          #                              .packages = c('data.table','survival'),
                          #                              .export = c('n.type.events','dt.tv', 'post.op.covid.model','newdata.pred')) %dopar% {
                          #                              cuminc.cox(n.type.events = n.type.events,
                          #                                         dt = 'dt.tv[patient_id %in% sample(unique(patient_id), replace = T) & (postop.covid.cohort)]', 
                          #                                         model = 'post.op.covid.model', 
                          #                                         newdata = 'newdata.pred',
                          #                                         day = 30)}             
                          #  t.samples <- t(apply(samples,1,quantile,c(0.25,0.5,0.75)))
                          #  boot.IQR <-apply(t.samples,1,function(x) paste0(x[2],' (',x[1],',',x[3],')'))

                           death.risk.30day <- predict(object = post.op.mort.model[[1]], 
                                                       newdata = newdata.pred, type = 'expected',se.fit = T)
                           

                           cbind(cuminc.cox(n.type.events = n.type.events,
                                                                  dt = 'dt.tv[start >=0 & any.op == T]', 
                                                                  model = 'post.op.mort.model', 
                                                                  newdata = 'newdata.pred',
                                                                  day = 30),

                           matrix(paste0(round((1- exp(-death.risk.30day$fit))*100,3),
                                         ' (', round((1 - exp(-(death.risk.30day$fit - 1.96*death.risk.30day$se.fit)))*100,3),',',
                                         round((1 - exp(-(death.risk.30day$fit + 1.96*death.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows)
                           )
                         })


save(post.op.mort.model,adjusted.mort.cuminc, file = here::here("output","postopmort_adjusted.RData"))

data.table::fwrite(adjusted.mort.cuminc, file = here::here("output","postopmort_adjusted.csv"))
