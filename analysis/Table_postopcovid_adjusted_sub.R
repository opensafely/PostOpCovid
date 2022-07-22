library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures.sub <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')
covariates <- c(procedures.sub,'sex','age.cat','bmi.cat','imd5','wave','LOS.bin',
                'vaccination.status.factor','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID','region')

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','postop.covid.cohort','end.fu'))]

dt.tv[,(drop.vars) := NULL]

dt.tv[,(procedures.sub) := lapply(.SD,function(x) x==1), .SDcols = (procedures.sub)]



data.table::setkey(dt.tv,patient_id,tstart,tstop)

n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
      (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))

post.op.covid.model.waves.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Colectomy + Cholecystectomy + KneeReplacement + age.cat + sex + bmi.cat + imd5 + 
                                                      vaccination.status.factor + region + Current.Cancer + Emergency*wave + LOS.bin + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort) & sub.op == T], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.waves.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelwavessub.csv"))


new.data.postop <- data.table::data.table(
  'start' = rep(0,8*length(procedures.sub)),
  'end' = rep(30,8*length(procedures.sub)),
  'event' = rep(F,8*length(procedures.sub)),
  'Colectomy' = c(rep(T,8),rep(F,24)),
  'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
  'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
  'KneeReplacement'=c(rep(F,24),rep(T,8)),
  'age.cat' = rep('(50,70]',8*length(procedures.sub)),
  'sex' = rep('F',8*length(procedures.sub)),
  'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures.sub)),
  'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures.sub)),
  'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures.sub)),
  'vaccination.status.factor' = rep('3',8*length(procedures.sub)),
  'region' = rep("East Midlands",8*length(procedures.sub)),
  'Current.Cancer' = rep(T,8*length(procedures.sub)),
  'Emergency' =  rep(c(rep(F,4),rep(T,4)), times = length(procedures.sub)),
  'LOS.bin' = rep(F,8*length(procedures.sub)),
  'Charl12' =  rep('Single',8*length(procedures.sub)),
  'recentCOVID' = rep(F,8*length(procedures.sub)),
  'previousCOVID' = rep(F,8*length(procedures.sub)),
  'patient_id' = 1:(8*length(procedures.sub)))

cuminc.adjusted.waves.sub <- 
  matrix(cuminc.cox(n.type.events = n.type.events,
                    dt = 'dt.tv[(postop.covid.cohort) & sub.op == T]',
                    model = 'post.op.covid.model.waves.sub', 
                    newdata = 'new.data.postop',
                    day = 30), byrow = T, ncol = 4)

colnames(cuminc.adjusted.waves.sub) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.waves.sub) <- paste0(c('Elective','Emergency'),rep(procedures.sub, each = 2))

data.table::fwrite(cuminc.adjusted.waves.sub, file = here::here("output","postopcovid_adjusted_waves_sub.csv"))


adjusted.waves.sub.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.waves.sub, keep.rownames = T),
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

ggplot2::ggsave(plot = adjusted.waves.sub.plot, here::here('output','adjusted_waves_sub_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

#############################################################################################

data.table::setkey(dt.tv,"patient_id","tstart","tstop")
post.op.covid.model.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~ Colectomy + Cholecystectomy  + KneeReplacement +
                                                      age.cat + sex + bmi.cat + imd5 + 
                                                      vaccination.status.factor + region + Current.Cancer + 
                                                      Emergency + wave + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                                    id = patient_id,
                                                    data = dt.tv[(postop.covid.cohort)  & sub.op == T], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelsub.csv"))


adjusted.cuminc.sub <-  data.table::as.data.table(foreach::foreach(predi = 1:length(covariates), .combine = 'rbind', .inorder = T) %do% {
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
                            if ( predi <= length(procedures.sub)) {
                                newdata.pred[,(procedures.sub) := F]
                                newdata.pred[,(procedures.sub[predi]) := c(F,T)]
                            } else {

                            
                          #  newdata.pred <- data.table::data.table('start' = rep(0,newdata.rows),
                          #                                         'end' = rep(30,newdata.rows),
                          #                                         'event' = rep(F,newdata.rows),
                          #                                         'patient_id' = 1:newdata.rows)
                          #  if ( predi > length(procedures.sub)) {
                          #    newdata.pred[,(procedures.sub) := lapply(procedures.sub, function(x) x == procedures.sub[which.max(dt.tv[,lapply(.SD,sum,na.rm = T), .SDcols = c(procedures.sub)])])] } else {
                          #      newdata.pred[,(procedures.sub) := lapply(procedures.sub, function(x) x == covariates[predi] & patient_id > 1)]
                          #    }
                           
                          #  newdata.pred[,(covariates[-c(1:length(procedures.sub))]) := lapply(((length(procedures.sub)+1):length(covariates)), function(i.c) {
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

                           
                           
                           
                           death.risk.30day <- predict(object = post.op.covid.model.sub[[3]], 
                                                       newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           readmit.risk.30day <- predict(object = post.op.covid.model.sub[[2]], 
                                                         newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           covid.risk.30day <- predict(object = post.op.covid.model.sub[[1]], 
                                                       newdata = newdata.pred,, type = 'expected',se.fit = T)
                           
                           cbind(matrix(paste0(round((1- exp(-covid.risk.30day$fit))*100,3),
                                               ' (', round((1 - exp(-(covid.risk.30day$fit - 1.96*covid.risk.30day$se.fit)))*100,3),',',
                                               round((1 - exp(-(covid.risk.30day$fit + 1.96*covid.risk.30day$se.fit)))*100,3),')'),nrow =newdata.rows),
                                 cuminc.cox(n.type.events = n.type.events,
                                            dt = 'dt.tv[(postop.covid.cohort) & sub.op == T]', 
                                            model = 'post.op.covid.model.sub', 
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


save(cuminc.adjusted.waves.sub,post.op.covid.model.sub,adjusted.cuminc.sub, file = here::here("output","postopcovid_adjusted_sub.RData"))
# Take out baseline no procedures.sub groups that are not observeds
data.table::fwrite(adjusted.cuminc.sub, file = here::here("output","postopcovid_adjusted_sub.csv"))
