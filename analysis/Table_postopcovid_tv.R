library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv.splits <- data.table::setDT(arrow::read_feather(here::here("output","cohort_postdisch_week_splits.feather")))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

n.type.events <- 1:2 #sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)

dt.tv.splits[event == 3, event := 2]
data.table::setkey(dt.tv.splits, patient_id, end.fu, start)
post.op.covid.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~  Abdominal+ Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + 
  sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency*week.post.disch  + Charl12 + recentCOVID + previousCOVID,
   id = patient_id,
  data = dt.tv.splits[(postop.covid.cohort) & start <=end], model = T))

data.table::fwrite(broom::tidy(post.op.covid.model.split[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelsplit.csv"))


newdata.rows <- length(levels(dt.tv.splits$week.post.disch)) - 1

newdata.pred <- data.table::data.table('start' = c(-7,0,7,14,21),
                                       'end' = c(0,7,14,21,28),
                                       'event' = rep(F,newdata.rows),
                                      'week.post.disch' = paste(0:(newdata.rows - 1)),
                                      'patient_id' = 1:newdata.rows,
                                      'Abdominal' = rep(T,newdata.rows),
                                      'Cardiac'= rep(F,newdata.rows),
                                      'Obstetrics'=rep(F,newdata.rows),
                                      'Orthopaedic'=rep(F,newdata.rows),
                                      'Thoracic'=rep(F,newdata.rows),
                                      'Vascular'=rep(F,newdata.rows),
                                      'age.cat' = rep('(50,70]',newdata.rows),
                                      'sex' = rep('F',newdata.rows),
                                      'bmi.cat' = rep(levels(dt.tv.splits$bmi.cat)[2],newdata.rows),
                                      'imd5' = rep(levels(dt.tv.splits$imd5)[3], newdata.rows),
                                      'wave' = rep(paste0('Wave_',4),times = newdata.rows),
                                      'vaccination.status.factor' = rep('3',newdata.rows),
                                      'region' = rep("East Midlands",newdata.rows),
                                      'Current.Cancer' = rep(T,newdata.rows),
                                      'Emergency' =  rep(F,newdata.rows),
                                      'Charl12' =  rep('Single',newdata.rows),
                                      'recentCOVID' = rep(F,newdata.rows),
                                      'previousCOVID' = rep(F,newdata.rows)
                                      )
#newdata.pred[,(procedures) := lapply(procedures, function(x) x == procedures[which.max(dt.tv[,lapply(.SD,sum,na.rm = T), .SDcols = c(procedures)])])]

# newdata.pred[,(covariates[-c(1:length(procedures))]) := lapply(((length(procedures)+1):length(covariates)), function(i.c) {
#   if(is.factor(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
#     as.character(rep(max.category(i.c),newdata.rows))
#   } else if(is.logical(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
#     as.logical(rep(max.category(i.c),newdata.rows))
#   } else if(is.numeric(dt.tv[!is.na(get(covariates[i.c])),get(covariates[i.c])])) {
#     is.numeric(rep(max.category(i.c),newdata.rows))
#   } else {
#     rep(max.category(i.c),newdata.rows)
#   }
# })]

base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.covid.model.split[[i]],centered = F))

base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
                                                                            base.haz = base.haz[[i]][,1] - 
                                                                              c(0,head(base.haz[[i]][,1],-1)))})

lp <- lapply(n.type.events, function(i) {  data.table::data.table('time' = seq(-7,21,7),
                                                                    'risk' = exp(predict(object = post.op.covid.model.split[[i]],
                                                                                         type = 'lp', 
                                                                                         newdata = newdata.pred))) })

base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))

for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[[j]])), j, 0)

weekly.post.op.risk <- 
  unlist(round(100*apply(exp(apply(safelog(1 - Reduce('+',lapply(n.type.events, function(i) {
    lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] 
  }))),2,cumsum))*
    lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] ,2,cumsum), digits = 3))

weekly.post.op.risk[!is.finite(weekly.post.op.risk)] <- 0

times.comb <- unique(sort(unlist(base.haz.merge$time)))[unique(sort(unlist(base.haz.merge$time))) >= -7]

weekly.post.op.risk <- c(weekly.post.op.risk[max(which(times.comb <= 0))],
                        weekly.post.op.risk[max(which(times.comb <= 7))],
                        weekly.post.op.risk[max(which(times.comb <= 14))],
                        weekly.post.op.risk[max(which(times.comb <= 21))],
                        weekly.post.op.risk[max(which(times.comb <= 28))])

weekly.post.op.risk[!is.finite(weekly.post.op.risk)] <- 0

weekly.post.op.risk  <-  data.table::data.table("Risk" = weekly.post.op.risk - c(0,weekly.post.op.risk[-length(weekly.post.op.risk)]),
                                               "Risk period" = c("Week pre discharge","1st week","2nd week","3rd week","4th week"))

# weekly.post.op.risk  <-  data.table::data.table("Risk" = (weekly.post.op.risk - c(0,weekly.post.op.risk[-length(weekly.post.op.risk)])[which(times.comb <= 28)]),
#                                                 "Days.post.discharge" = (-7):28)

data.table::fwrite(weekly.post.op.risk, file = here::here("output","postopcovid_tv.csv"))

##############
post.op.VTE.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event.VTE==i) ~  Abdominal + survival::strata(postcovid) +    
                                            Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + 
                                            sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + 
                                            Current.Cancer + Emergency*week.post.disch + Charl12 + recentCOVID + previousCOVID, 
                                          id = patient_id,
                                          data = dt.tv.splits[(postcovid.VTE.cohort) & start <=end], model = T))

data.table::fwrite(broom::tidy(post.op.VTE.model.split[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopVTEmodelsplit.csv"))


newdata.rows <- length(levels(dt.tv.splits$week.post.disch)) - 1

newdata.pred <- data.table::data.table('start' = rep(c(-7,0,7,14,21), times = 2),
                                       'end' = rep(c(0,7,14,21,28),times = 2),
                                       'event.VTE' = rep(F,newdata.rows*2),
                                       'week.post.disch' = rep(paste(0:(newdata.rows - 1)), times = 2),
                                       'patient_id' = rep(1:2,each = newdata.rows),
                                       'Abdominal' = rep(T,newdata.rows*2),
                                       'Cardiac'= rep(F,newdata.rows*2),
                                       'Obstetrics'=rep(F,newdata.rows*2),
                                       'Orthopaedic'=rep(F,newdata.rows*2),
                                       'Thoracic'=rep(F,newdata.rows*2),
                                       'postcovid'=rep(c(0,1), each = newdata.rows),
                                       'Vascular'=rep(F,newdata.rows*2),
                                       'age.cat' = rep('(50,70]',newdata.rows*2),
                                       'sex' = rep('F',newdata.rows*2),
                                       'bmi.cat' = rep(levels(dt.tv.splits$bmi.cat)[2],newdata.rows*2),
                                       'imd5' = rep(levels(dt.tv.splits$imd5)[3], newdata.rows*2),
                                       'wave' = rep(paste0('Wave_',4),times = newdata.rows*2),
                                       'vaccination.status.factor' = rep('3',newdata.rows*2),
                                       'region' = rep("East Midlands",newdata.rows*2),
                                       'Current.Cancer' = rep(T,newdata.rows*2),
                                       'Emergency' =  rep(F,newdata.rows*2),
                                       'Charl12' =  rep('Single',newdata.rows*2),
                                       'recentCOVID' = rep(F,newdata.rows*2),
                                       'previousCOVID' = rep(F,newdata.rows*2)
)

base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.VTE.model.split[[i]],centered = F))

base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
                                                                            base.haz = base.haz[[i]][,1] - 
                                                                              c(0,head(base.haz[[i]][,1],-1)))})


lp <- lapply(n.type.events, function(i) {  data.table::dcast(data.table::data.table('patient_id' = rep(1:2, each = 5),
  'time' = rep(seq(-7,21,7),2),
                                                                    'risk' = exp(predict(object = post.op.VTE.model.split[[i]],
                                                                                         type = 'lp', 
                                                                                         newdata = newdata.pred))),time ~patient_id, value.var = 'risk')})

base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))

for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[[j]])), j, 0)

weekly.post.op.VTE.risk <- 
  unlist(round(100*apply(exp(apply(safelog(1 - Reduce('+',lapply(n.type.events, function(i) {
    lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)] 
  }))),2,cumsum))*
    lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)], 2, cumsum ), digits = 3))

for (j in 1:ncol(weekly.post.op.VTE.risk)) set(weekly.post.op.VTE.risk, which(!is.finite(weekly.post.op.VTE.risk[[j]])), j, 0)

times.comb <- unique(sort(unlist(base.haz.merge$time)))[unique(sort(unlist(base.haz.merge$time))) >= -7]

weekly.post.op.VTE.risk <- rbind(weekly.post.op.VTE.risk[max(which(times.comb <= 0)),],
                            weekly.post.op.VTE.risk[max(which(times.comb <= 7)),],
                            weekly.post.op.VTE.risk[max(which(times.comb <= 14)),],
                            weekly.post.op.VTE.risk[max(which(times.comb <= 21)),],
                        weekly.post.op.VTE.risk[max(which(times.comb <= 28)),])

for (j in 1:ncol(weekly.post.op.VTE.risk)) set(weekly.post.op.VTE.risk, which(!is.finite(weekly.post.op.VTE.risk[[j]])), j, 0)

weekly.post.op.VTE.risk <- as.vector(weekly.post.op.VTE.risk[,lapply(.SD, function(x) as.numeric(x) - as.numeric(c(0,head(x,-1))))])
weekly.post.op.VTE.risk  <-  cbind(data.table::data.table("COVID"= rep(c(F,T), each = 5),weekly.post.op.VTE.risk,c("Week pre discharge","1st week","2nd week","3rd week","4th week")))


##################################
save(weekly.post.op.risk,weekly.post.op.VTE.risk, file = here::here("output","postopcovid_tv.RData"))
data.table::fwrite(weekly.post.op.VTE.risk, file = here::here("output","postopcovid_VTE_tv.csv"))


