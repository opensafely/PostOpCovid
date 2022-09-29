library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

#dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_postdisch_week_splits.feather")))
#dt.tv <- dt.tv[,.(patient_id,Abdominal,Cardiac,Obstetrics,Thoracic,Vascular,age.cat,sex,bmi.cat,imd5,wave,vaccination.status.factor,region,Current.Cancer,Emergency,
#                                week.post.disch,week.post.op,LOS.bin,Charl12,recentCOVID,previousCOVID,postcovid,start.readmit,end.readmit,
#                                tstart,tstop,end.fu,start,end,event,postop.covid.cohort,los.end,event.VTE,postcovid.VTE.cohort,study.start)]
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

# Start = 0 = day of operation
#dt.tv[, `:=`(start = tstart - study.start,
#                    end = tstop - study.start)]
#dt.tv <- dt.tv[start >= 0,] # Need to start follow up on day after operation as can't identify order when events on same day

# Not enough deaths to treat separately from emergency readmissions
#dt.tv[event == 3, event := 2]
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

data.table::setkey(dt.tv,patient_id,tstart,tstop)


data.table::setkey(dt.tv, patient_id, end.fu, start)
post.op.covid.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~  Abdominal  + Obstetrics +  CardioThoracicVascular + age.cat + 
  sex  + imd5 + wave + vaccination.status.factor  + Current.Cancer + Emergency +  Charl12 + recentCOVID + previousCOVID,
   id = patient_id,
  data = dt.tv[(postop.covid.cohort) & start <=end], model = T))

#data.table::fwrite(broom::tidy(post.op.covid.model.split[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopcovidmodelsplit.csv"))
names(post.op.covid.model.split) <- c('COVID-19','Non COVID-19 emergency readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.covid.model.split,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, conf_level = .95, exponentiate = TRUE, output = here::here("output","postopcovidmodelsplit.html"))
newdata.rows <- 1

newdata.pred <- data.table::data.table('start' = c(0),
                                       'end' = c(35),
                                       'event' = rep(F,newdata.rows),
                                      'week.post.op' = paste(1:newdata.rows),
                                      'patient_id' = 1:newdata.rows,
                                      'Abdominal' = rep(T,newdata.rows),
                                    #  'Cardiac'= rep(F,newdata.rows),
                                      'Obstetrics'=rep(F,newdata.rows),
                                      'Orthopaedic'=rep(F,newdata.rows),
                                     # 'Thoracic'=rep(F,newdata.rows),
                                      'CardioThoracicVascular'=rep(F,newdata.rows),
                                      'age.cat' = rep('(50,70]',newdata.rows),
                                      'sex' = rep('F',newdata.rows),
                                      'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],newdata.rows),
                                      'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows),
                                      'wave' = rep(paste0('Wave_',3),times = newdata.rows),
                                      'vaccination.status.factor' = rep('3',newdata.rows),
                                      'region' = rep("East Midlands",newdata.rows),
                                      'Current.Cancer' = rep(T,newdata.rows),
                                      'Emergency' =  rep(F,newdata.rows),
                                      'LOS.bin' =  rep(F,newdata.rows),
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

# base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.covid.model.split[[i]],centered = F))
# 
# base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
#                                                                             base.haz = base.haz[[i]][,1] - 
#                                                                               c(0,head(base.haz[[i]][,1],-1)))})
# 
# lp <- lapply(n.type.events, function(i) {  data.table::data.table('time' = seq(0,28,7),#0:34,
#                                                                     'risk' = exp(predict(object = post.op.covid.model.split[[i]],
#                                                                                          type = 'lp', 
#                                                                                          newdata = newdata.pred))) })
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[[j]])), j, 0)
# 
# weekly.post.op.risk <- 
#   unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
#     lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] 
#   })),2,cumsum))*
#     lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] ,2,cumsum), digits = 3))
# 

weekly.post.op.risk <- cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postop.covid.cohort) & start <=end]', model = 'post.op.covid.model.split', newdata = 'newdata.pred', day = 1:35)

weekly.post.op.risk[!is.finite(weekly.post.op.risk)] <- 0

times.comb <-sort(unique(dt.tv[(postop.covid.cohort) & start <=end][ event %in% n.type.events ,end]))

#times.comb <- unique(sort(unlist(base.haz.merge$time)))[unique(sort(unlist(base.haz.merge$time))) >= 0]

weekly.post.op.risk <- c(weekly.post.op.risk[max(which(times.comb < 7))],
                        weekly.post.op.risk[max(which(times.comb < 14))],
                        weekly.post.op.risk[max(which(times.comb < 21))],
                        weekly.post.op.risk[max(which(times.comb < 28))],
                        weekly.post.op.risk[max(which(times.comb < 35))])

weekly.post.op.risk[!is.finite(weekly.post.op.risk)] <- 0

weekly.post.op.risk  <-  data.table::data.table("Risk" = weekly.post.op.risk - c(0,weekly.post.op.risk[-length(weekly.post.op.risk)]),
                                               "Risk period" = c("1st week","2nd week","3rd week","4th week","5th week"))

# weekly.post.op.risk  <-  data.table::data.table("Risk" = (weekly.post.op.risk - c(0,weekly.post.op.risk[-length(weekly.post.op.risk)])[which(times.comb <= 28)]),
#                                                 "Days.post.discharge" = (-7):28)

data.table::fwrite(weekly.post.op.risk, file = here::here("output","postopcovid_tv.csv"))

##############

# Not enough deaths to treat separately from emergency readmissions
#dt.tv[event.VTE == 3, event.VTE := 2]
n.type.events <- sort(unique(dt.tv[(postcovid.VTE.cohort),event.VTE]))[-1]


dt.tv <- dt.tv[start>0 & end.readmit <=90] # Need to start follow up day after discharge to avoid discharge diagnoses


post.op.VTE.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.VTE==i) ~  Abdominal + postcovid +    
                                         Obstetrics + CardioThoracicVascular+ age.cat + 
                                            sex + imd5 + wave + vaccination.status.factor  + 
                                            Current.Cancer + Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                          id = patient_id,
                                          data = dt.tv[(postcovid.VTE.cohort) & start <=end], model = T))

#data.table::fwrite(broom::tidy(post.op.VTE.model.split[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopVTEmodelsplit.csv"))
names(post.op.VTE.model.split) <- c('Post discharge VTE','Non COVID-19 emergency readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.VTE.model.split,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodelsplit.html"))


newdata.rows <- 1 #length(levels(dt.tv$week.post.disch)) - 1

newdata.pred <- data.table::data.table('start.readmit' = rep(c(-7), times = 2),
                                       'end.readmit' = rep(c(35),times = 2),
                                       'event.VTE' = rep(F,newdata.rows*2),
                                       'week.post.disch' = rep(paste(0:(newdata.rows - 1)), times = 2),
                                       'patient_id' = rep(1:2,each = newdata.rows),
                                       'Abdominal' = rep(T,newdata.rows*2),
                                   #    'Cardiac'= rep(F,newdata.rows*2),
                                       'Obstetrics'=rep(F,newdata.rows*2),
                                       'Orthopaedic'=rep(F,newdata.rows*2),
                                       'CardioThoracicVascular'=rep(F,newdata.rows*2),
                                    #   'Thoracic'=rep(F,newdata.rows*2),
                                       'postcovid'=rep(c(0,1), each = newdata.rows),
                                       'age.cat' = rep('(50,70]',newdata.rows*2),
                                       'sex' = rep('F',newdata.rows*2),
                                       'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],newdata.rows*2),
                                       'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows*2),
                                       'wave' = rep(paste0('Wave_',3),times = newdata.rows*2),
                                       'vaccination.status.factor' = rep('3',newdata.rows*2),
                                       'region' = rep("East Midlands",newdata.rows*2),
                                       'Current.Cancer' = rep(T,newdata.rows*2),
                                       'Emergency' =  rep(F,newdata.rows*2),
                                       'LOS.bin' =  rep(F,newdata.rows*2),
                                       'Charl12' =  rep('Single',newdata.rows*2),
                                       'recentCOVID' = rep(F,newdata.rows*2),
                                       'previousCOVID' = rep(F,newdata.rows*2)
)

# base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.VTE.model.split[[i]],centered = F))
# 
# base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
#                                                                             base.haz = base.haz[[i]][,1] - 
#                                                                               c(0,head(base.haz[[i]][,1],-1)))})
# 
# 
# lp <- lapply(n.type.events, function(i) {  data.table::dcast(data.table::data.table('patient_id' = rep(1:2, each = 5),
#   'time' = rep(seq(0,28,7),2),
#                                                                     'risk' = exp(predict(object = post.op.VTE.model.split[[i]],
#                                                                                          type = 'lp', 
#                                                                                          newdata = newdata.pred))),time ~patient_id, value.var = 'risk')})
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[,(j)])), j, 0)
# 

weekly.post.op.VTE.risk <- cuminc.cox(n.type.events = n.type.events,
                                      dt = 'dt.tv[(postcovid.VTE.cohort) & start <=end]', 
                                      model = 'post.op.VTE.model.split', newdata = 'newdata.pred', day = -7:35)


  # unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
  #   lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)] 
  # })),2,cumsum))*
  #   lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)], 2, cumsum ), digits = 3))

weekly.post.op.VTE.risk[!is.finite(weekly.post.op.VTE.risk)] <- 0

times.comb <-sort(unique(dt.tv[(postcovid.VTE.cohort) & start <=end][ event.VTE %in% n.type.events ,end]))


weekly.post.op.VTE.risk <- rbind(weekly.post.op.VTE.risk[max(which(times.comb <= 0)),],
                                 weekly.post.op.VTE.risk[max(which(times.comb <= 7)),],
                            weekly.post.op.VTE.risk[max(which(times.comb <= 14)),],
                            weekly.post.op.VTE.risk[max(which(times.comb <= 21)),],
                        weekly.post.op.VTE.risk[max(which(times.comb <= 28)),],
                        weekly.post.op.VTE.risk[max(which(times.comb <= 35)),])

weekly.post.op.VTE.risk[!is.finite(weekly.post.op.VTE.risk)] <- 0

weekly.post.op.VTE.risk <- as.vector(weekly.post.op.VTE.risk - rbind(c(0,0),weekly.post.op.VTE.risk[1:5,]))
weekly.post.op.VTE.risk  <-  cbind(data.table::data.table("COVID"= rep(c(F,T), each = 6),
"Risk" = weekly.post.op.VTE.risk,
"Risk period" = c("Week in hospital","1st post discharge week","2nd week","3rd week","4th week","5th week")))


##################################
save(weekly.post.op.risk,weekly.post.op.VTE.risk, file = here::here("output","postopcovid_tv.RData"))
data.table::fwrite(weekly.post.op.VTE.risk, file = here::here("output","postopcovid_VTE_tv.csv"))


