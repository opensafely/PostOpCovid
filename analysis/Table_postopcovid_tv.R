library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))
max.category <- function(predi) {unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)[which.max(dt.tv[!is.na(get(covariates[predi])),.N, by = c(covariates[predi])][['N']])]}

###########################################################

load(file = here::here("output","cohort_long.RData"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')


data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer',
                'Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv[,`:=`(post.disch.wk1 = discharge.date  ,
            post.disch.wk2 = discharge.date + 7 ,
            post.disch.wk3 = discharge.date + 14 ,
            post.disch.wk4 = discharge.date + 21 ,
            post.disch.wk5 = discharge.date + 28,
            covid.event = event == 1,
            censor.event = event!=1)]


dt.tv.splits <- data.table::melt(dt.tv[is.finite(discharge.date),tail(.SD,1),
                                       .SDcols = c(paste0("post.disch.wk",1:5)),
                                       keyby = c('patient_id')], 
                                 id.vars = c('patient_id'),
                                 value.name = "post.disch.wk",
                                 variable.name = "week.post.disch" )[
                                   order(patient_id, week.post.disch),][,
                                        `:=` (tstart = post.disch.wk,
                                              week.post.disch = as.numeric(gsub(".*?([0-9]+).*",
                                                                                '\\1',
                                                                                week.post.disch)))][
                                                ,.(patient_id, week.post.disch,tstart)]

#
dt.tv[,week.post.disch := NA]
dt.tv.splits <- unique(data.table::rbindlist(list(dt.tv.splits,dt.tv[,.(patient_id, week.post.disch, tstart)])))
dt.tv.splits[,post.disch.wk := tstart]
data.table::setkey(dt.tv.splits,patient_id, tstart)
dt.tv[,week.post.disch := NULL]
dt.tv.splits <- dt.tv[dt.tv.splits,,roll = Inf, on = .(patient_id,tstart)]

dates.expand.start.align_(dt = 'dt.tv.splits',
                          start.DTTM = 'tstart',
                          end.DTTM = 'tstop',
                          ID = 'patient_id',
                          merged.DTTM = 'post.disch.wk')

dates.expand.end.align_(dt = 'dt.tv.splits',
                        start.DTTM = 'tstart',
                        end.DTTM = 'tstop',
                        ID = 'patient_id',
                        merged.DTTM = 'post.disch.wk')

locf.roll_(dt = 'dt.tv.splits',
           ID = 'patient_id',
           start.DTTM = 'tstart',
           group = 'c("patient_id","end.fu")',
           var.cols = 'c("week.post.disch")')

dt.tv.splits[tstart >= admit.date & (!is.finite(discharge.date) | tstop <= discharge.date), week.post.disch := 0 ]

dt.tv.splits[,week.post.disch := as.factor(week.post.disch)]
dt.tv.splits[, los.end := min(los.end, na.rm = T), keyby = .(patient_id, end.fu)]

dt.tv.splits[, `:=`(start = tstart - los.end,
                    end = tstop - los.end)]
dt.tv.splits <- dt.tv.splits[is.finite(los.end)]
dt.tv.splits[event == 3, event := 2]
data.table::setkey(dt.tv.splits, patient_id, end.fu, start)
post.op.covid.model.split <- 
  lapply(1:2, function(i) survival::coxph(survival::Surv(start,end,event==i) ~  Abdominal*week.post.disch + Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,data = dt.tv.splits[(postop.covid.cohort) & start <=end], model = T))

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
n.type.events <- 2

times <- lapply(1:n.type.events, function(i) { survival::basehaz(post.op.covid.model.split[[i]],centered = F)$time })

base.haz <- lapply(1:n.type.events, function(i) { data.table::data.table('time' = times[[i]],
                                                                         'base.haz' = survival::basehaz(post.op.covid.model.split[[i]],centered = F)[,1] -
                                                                           c(0,head(survival::basehaz(post.op.covid.model.split[[i]],centered = F)[,1],-1)))})

lp <- lapply(1:n.type.events, function(i) {  data.table::data.table('time' = seq(-7,21,7),
                                                                    'risk' = exp(predict(object = post.op.covid.model.split[[i]],
                                                                                         type = 'lp', 
                                                                                         newdata = newdata.pred))) })

base.haz.merge <- Reduce(x =base.haz,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))


weekly.post.op.risk <- 
  unlist(round(100*cumsum(exp(cumsum(safelog(1 - Reduce('+',lapply(n.type.events, function(i) {
    lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.SD[,1]*.SD[,2] ,.SDcols = c(2:3)] 
  })))))*
    lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.SD[,1]*.SD[,2] ,.SDcols = c(2:3)] ), digits = 3))

weekly.post.op.risk[!is.finite(weekly.post.op.risk)] <- 0

times.comb <- unique(sort(unlist(times)))[unique(sort(unlist(times))) >= -7]

weekly.post.op.risk <- c(weekly.post.op.risk[max(which(times.comb <= 0))],
                         weekly.post.op.risk[max(which(times.comb <= 7))],
                         weekly.post.op.risk[max(which(times.comb <= 14))],
                         weekly.post.op.risk[max(which(times.comb <= 21))],
                         weekly.post.op.risk[max(which(times.comb <= 28))])


weekly.post.op.risk  <-  data.table::data.table("Risk" = weekly.post.op.risk - c(0,weekly.post.op.risk[-length(weekly.post.op.risk)]),
                                                "Risk period" = c("Week pre discharge","1st week","2nd week","3rd week","4th week"))

##############
post.op.VTE.model.split <- 
  lapply(1:2, function(i) survival::coxph(survival::Surv(start,end,event.VTE==i) ~  Abdominal + survival::strata(postcovid)*week.post.disch +  
                                            Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + age.cat + 
                                            sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + 
                                            Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, 
                                          id = patient_id,
                                          data = dt.tv.splits[(postcovid.VTE.cohort) & start <=end], model = T))

data.table::fwrite(broom::tidy(post.op.VTE.model.split[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopVTEmodelsplit.csv"))


newdata.rows <- length(levels(dt.tv.splits$week.post.disch)) - 1

newdata.pred <- data.table::data.table('start' = rep(c(-7,0,7,14,21), times = 2),
                                       'end' = rep(c(0,7,14,21,28),times = 2),
                                       'event' = rep(F,newdata.rows*2),
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
                                       'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],newdata.rows*2),
                                       'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows*2),
                                       'wave' = rep(paste0('Wave_',4),times = newdata.rows*2),
                                       'vaccination.status.factor' = rep('3',newdata.rows*2),
                                       'region' = rep("East Midlands",newdata.rows*2),
                                       'Current.Cancer' = rep(T,newdata.rows*2),
                                       'Emergency' =  rep(F,newdata.rows*2),
                                       'Charl12' =  rep('Single',newdata.rows*2),
                                       'recentCOVID' = rep(F,newdata.rows*2),
                                       'previousCOVID' = rep(F,newdata.rows*2)
)

n.type.events <- 2

times <- lapply(1:n.type.events, function(i) { survival::basehaz(post.op.VTE.model.split[[i]],centered = F)$time })

base.haz <- lapply(1:n.type.events, function(i) { data.table::data.table('time' = times[[i]],
                                                                         'base.haz' = survival::basehaz(post.op.VTE.model.split[[i]],centered = F)[,1] -
                                                                           c(0,head(survival::basehaz(post.op.VTE.model.split[[i]],centered = F)[,1],-1)))})

lp <- lapply(1:n.type.events, function(i) {  data.table::dcast(data.table::data.table('patient_id' = rep(1:2, each = 5),
  'time' = rep(seq(-7,21,7),2),
                                                                    'risk' = exp(predict(object = post.op.VTE.model.split[[i]],
                                                                                         type = 'lp', 
                                                                                         newdata = newdata.pred))),time ~patient_id, value.var = 'risk')})

base.haz.merge <- Reduce(x =base.haz,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))


weekly.post.op.VTE.risk <- 
  unlist(round(100*apply(exp(apply(safelog(1 - Reduce('+',lapply(n.type.events, function(i) {
    lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)] 
  }))),2,cumsum))*
    lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= -7][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)], 2, cumsum ), digits = 3))

weekly.post.op.VTE.risk[!is.finite(weekly.post.op.VTE.risk)] <- 0

times.comb <- unique(sort(unlist(times)))[unique(sort(unlist(times))) >= -7]

weekly.post.op.VTE.risk <- rbind(weekly.post.op.VTE.risk[max(which(times.comb <= 0)),],
                             weekly.post.op.VTE.risk[max(which(times.comb <= 7)),],
                             weekly.post.op.VTE.risk[max(which(times.comb <= 14)),],
                             weekly.post.op.VTE.risk[max(which(times.comb <= 21)),],
                         weekly.post.op.VTE.risk[max(which(times.comb <= 28)),])


weekly.post.op.VTE.risk  <-  data.table::data.table("COVID"= rep(c(F,T), each = 5),
  "Risk" = as.vector(weekly.post.op.VTE.risk - rbind(c(0,0),weekly.post.op.VTE.risk[-nrow(weekly.post.op.VTE.risk),])),
                                                "Risk period" = c("Week pre discharge","1st week","2nd week","3rd week","4th week"))

##################################
save(weekly.post.op.risk,weekly.post.op.VTE.risk, file = here::here("output","postopcovid_tv.RData"))
