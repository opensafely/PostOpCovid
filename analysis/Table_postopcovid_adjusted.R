library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')


covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','postop.covid.cohort','end.fu'))]

dt.tv[,(drop.vars) := NULL]
dt.tv <- dt.tv[(postop.covid.cohort)]
gc()
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

post.op.covid.model.waves <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Abdominal*wave +  Obstetrics*wave + CardioThoracicVascular*wave + age.cat + sex + bmi.cat + imd5 + region + 
                                                      vaccination.status.factor + Current.Cancer + Emergency*wave + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort)], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.waves[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelwaves.csv"))


new.data.postop <- data.table::data.table(
  'start' = rep(0,8*length(procedures)),
  'end' = rep(30,8*length(procedures)),
  'event' = rep(F,8*length(procedures)),
  'Abdominal' = c(rep(T,8),rep(F,24)),
 # 'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
  'Obstetrics'=c(rep(F,8),rep(T,8),rep(F,16)),
  'Orthopaedic'=c(rep(F,16),rep(T,8),rep(F,8)),
  #'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
  'CardioThoracicVascular'=c(rep(F,24),rep(T,8)),
  'age.cat' = rep('(50,70]',8*length(procedures)),
  'sex' = rep('F',8*length(procedures)),
  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
  'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
  'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
  'vaccination.status.factor' = rep('3',8*length(procedures)),
  'region' = rep("East Midlands",8*length(procedures)),
  'Current.Cancer' = rep(T,8*length(procedures)),
  'Emergency' =  rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
  'LOS.bin' = rep(F,8*length(procedures)),
  'Charl12' =  rep('Single',8*length(procedures)),
  'recentCOVID' = rep(F,8*length(procedures)),
  'previousCOVID' = rep(F,8*length(procedures)),
  'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.waves <- 
  matrix(cuminc.cox(n.type.events = n.type.events,
                    dt = 'dt.tv[(postop.covid.cohort)]',
                    model = 'post.op.covid.model.waves', 
                    newdata = 'new.data.postop',
                    day = 30), byrow = T, ncol = 4)

colnames(cuminc.adjusted.waves) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.waves) <- paste0(c('Elective','Emergency'),rep(procedures, each = 2))

data.table::fwrite(cuminc.adjusted.waves, file = here::here("output","postopcovid_adjusted_waves.csv"))



adjusted.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.waves, keep.rownames = T),
                                 id.vars = 'rn',
                                 variable.name = 'Wave',
                                 value.name = '30 Day COVID Cumulative Incidence (%)')[, `:=`(Emergency = grepl('Emergency*',rn),
                                                                                          Operation = gsub('Emergency|Elective', '',rn))],
                ggplot2::aes(x = Wave, 
                             y = `30 Day COVID Cumulative Incidence (%)`, 
                             group = rn,
                             colour = Operation,
                             linetype = Emergency)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = adjusted.waves.plot, here::here('output','adjusted_waves_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

#############################################################################################

data.table::setkey(dt.tv,"patient_id","tstart","tstop")
post.op.covid.model <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Abdominal  + 
                                                      Obstetrics + CardioThoracicVascular + 
                                                      age.cat + sex  + bmi.cat + imd5 +  wave +  
                                                      vaccination.status.factor  + region +  Current.Cancer + 
                                                      Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                                    id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort)], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodel.csv"))


#library(doParallel)
#ncores <- parallel::detectCores(logical = F)
#cl <- parallel::makeCluster(ncores)
#doParallel::registerDoParallel(cl)

adjusted.cuminc <-  data.table::as.data.table(foreach::foreach(predi = 1:length(covariates),
                                                               .combine = 'rbind', 
                                                               .inorder = T) %do% {
                           newdata.rows <- length(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])]))
                           
   
                           newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                                                                  'end' = rep(30,newdata.rows),
                                                                  'event' = rep(F,newdata.rows),
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
                                                                  'wave' = rep(paste0('Wave_',3),times = newdata.rows),
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
                              newdata.pred[,(covariates[predi]) :=  as.character(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.logical(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  as.logical(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else if(is.numeric(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])])) {
                              newdata.pred[,(covariates[predi]) :=  is.numeric(sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)))]
                            } else {
                              newdata.pred[,(covariates[predi]) := sort(unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T))]
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

                           death.risk.30day <- predict(object = post.op.covid.model[[3]], 
                                                       newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           readmit.risk.30day <- predict(object = post.op.covid.model[[2]], 
                                                       newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           covid.risk.30day <- predict(object = post.op.covid.model[[1]], 
                                   newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           cbind(matrix(paste0(round((1- exp(-covid.risk.30day$fit))*100,3),
                                               ' (', round((1 - exp(-(covid.risk.30day$fit - 1.96*covid.risk.30day$se.fit)))*100,3),',',
                                               round((1 - exp(-(covid.risk.30day$fit + 1.96*covid.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows),
                                  cuminc.cox(n.type.events = n.type.events,
                                                                  dt = 'dt.tv[(postop.covid.cohort)]', 
                                                                  model = 'post.op.covid.model', 
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


save(cuminc.adjusted.waves,post.op.covid.model,adjusted.cuminc, file = here::here("output","postopcovid_adjusted.RData"))

data.table::fwrite(adjusted.cuminc, file = here::here("output","postopcovid_adjusted.csv"))
