library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

#dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_postdisch_week_splits.feather")))
#dt.tv <- dt.tv[,.(patient_id,Abdominal,Cardiac,Obstetrics,Thoracic,Vascular,age.cat,sex,bmi.cat,imd5,wave,vaccination.status.factor,region,Current.Cancer,Emergency,
#                                week.post.disch,week.post.op,LOS.bin,Charl12,recentCOVID,previousCOVID,postcovid,start.readmit, end.readmit,
#                                tstart,tstop,end.fu,start,end,event,postop.covid.cohort,los.end,event.VTE,postcovid.VTE.cohort,study.start)]
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
days.limit <- 35
# Start = 0 = day of operation
#dt.tv[, `:=`(start = tstart - study.start,
#                    end = tstop - study.start)]
#dt.tv <- dt.tv[start >= 0,] # Need to start follow up on day after operation as can't identify order when events on same day

#dt.tv[event == 3, event := 2]
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

data.table::setkey(dt.tv,patient_id,tstart,tstop)


data.table::setkey(dt.tv, patient_id, end.fu, start)
post.op.covid.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~  Abdominal + Cardiac + Obstetrics +  Thoracic + Vascular + age.cat + 
  sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency +  Charl12 + recentCOVID + previousCOVID,
   id = patient_id,
  data = dt.tv[(postop.covid.cohort) & start <=end], model = T))

newdata.rows <- 1 #

newdata.pred <- data.table::data.table('start'  = 0,
                                       'end' = days.limit,
                                       'event' = rep(F,newdata.rows),
                                  #    'week.post.op' = paste(rep(1:5,each = 7)),
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
# lp <- lapply(n.type.events, function(i) {  data.table::data.table('time' = 0:(newdata.rows - 1),
#                                                         'risk' = exp(predict(object = post.op.covid.model.split[[i]],
#                                                                              type = 'lp', 
#                                                                              newdata = newdata.pred))) })
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[[j]])), j, 0)
# 
# daily.post.op.risk <- 
#   unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
#     lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] 
#   })),2,cumsum))*
#     lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] ,2,cumsum), digits = 3))
daily.post.op.risk <- cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postop.covid.cohort) & start <=end]', model = 'post.op.covid.model.split', newdata = 'newdata.pred', day = 1:days.limit)

times.comb <-sort(unique(dt.tv[(postop.covid.cohort) & start <=end][ event %in% n.type.events & end <= days.limit ,end]))

daily.post.op.risk[!is.finite(daily.post.op.risk)] <- 0

#times.comb <- unique(sort(unlist(base.haz.merge$time)))[unique(sort(unlist(base.haz.merge$time))) >= 0]

daily.post.op.risk  <-  data.table::data.table(`risk` = daily.post.op.risk - c(0,daily.post.op.risk[-length(daily.post.op.risk)]),
                                                `Days post op` = times.comb)

data.table::fwrite(daily.post.op.risk, file = here::here("output","daily_postopcovid_tv.csv"))

ggplot2::ggplot(daily.post.op.risk[`Days post op`<= days.limit]) +
  ggplot2::geom_line(ggplot2::aes(x = `Days post op`, y = risk)) + 
  ggplot2::geom_smooth(ggplot2::aes(x = `Days post op`, y = risk)) +
  ggplot2::ylab("Daily risk (%)") +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle("Daily risk of COVID from day of operation", subtitle = "Predicted for abdominal elective cancer operation in 50-70 year old female\n BMI 20-25, 3rd IMD quintile with single morbidity\n booster vaccination and no previous COVD")
ggplot2::ggsave( filename = here::here("output","dailyCovidRisk.pdf"),device = "pdf",width = 8, height = 8, units = 'in',dpi = 'retina')
##############
# Not enough deaths to treat separately from emergency readmissions
#dt.tv[event.VTE == 3, event.VTE := 2]
n.type.events <- sort(unique(dt.tv[(postcovid.VTE.cohort),event.VTE]))[-1]

post.op.VTE.model.split <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.VTE==i) ~  Abdominal + postcovid +    
                                            Cardiac + Obstetrics + Thoracic + Vascular + age.cat + 
                                            sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + 
                                            Current.Cancer + Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                          id = patient_id,
                                          data = dt.tv[(postcovid.VTE.cohort) & start <=end], model = T))

newdata.rows = 1

newdata.pred <- data.table::data.table('start.readmit' = rep(0, times = 2),
                                       'end.readmit' = rep(days.limit,times = 2),
                                       'event.VTE' = rep(F,newdata.rows*2),
                                  #     'week.post.disch' = paste(rep(rep(1:5, each = 7), times = 2)),
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
# lp <- lapply(n.type.events, function(i) {  data.table::dcast(
#   data.table::data.table('patient_id' = rep(1:2, each = newdata.rows),
#   'time' = rep(0:(newdata.rows - 1),times = 2),
#               'risk' = exp(predict(object = post.op.VTE.model.split[[i]],
#                                    type = 'lp', 
#                                    newdata = newdata.pred))),time ~patient_id, value.var = 'risk')})
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[,(j)])), j, 0)

# daily.post.op.VTE.risk <- 
#   unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
#     lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)] 
#   })),2,cumsum))*
#     lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)], 2, cumsum ), digits = 3))

daily.post.op.VTE.risk <- cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & start <=end]', model = 'post.op.VTE.model.split', newdata = 'newdata.pred', day = 1:days.limit)

daily.post.op.VTE.risk[!is.finite(daily.post.op.VTE.risk)] <- 0

times.comb <-sort(unique(dt.tv[(postcovid.VTE.cohort) & start <=end][ event.VTE %in% n.type.events & end <=days.limit,end]))

daily.post.op.VTE.risk <- as.vector(daily.post.op.VTE.risk - rbind(c(0,0),daily.post.op.VTE.risk[-nrow(daily.post.op.VTE.risk),]))
daily.post.op.VTE.risk  <-  cbind(data.table::data.table("COVID"= rep(c(F,T), each = length(times.comb)),
"Risk" = daily.post.op.VTE.risk,
"Days post discharge" = rep(times.comb, times = 2)))

##################################
save(daily.post.op.risk,daily.post.op.VTE.risk, file = here::here("output","postopcovid_tv_daily.RData"))
data.table::fwrite(daily.post.op.VTE.risk, file = here::here("output","daily_postopcovid_VTE_tv.csv"))


######################
ggplot2::ggplot(daily.post.op.VTE.risk[`Days post discharge`<= newdata.rows]) +
  ggplot2::geom_line(ggplot2::aes(x = `Days post discharge`, y = Risk, colour = COVID)) + 
  ggplot2::geom_smooth(ggplot2::aes(x = `Days post discharge`, y = Risk, colour = COVID)) + 
  ggplot2::ylab("Daily risk (%)") +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle("Daily risk of VTE from day of discharge",
   subtitle = "Predicted for abdominal elective cancer operation in 50-70 year old female\n BMI 20-25, 3rd IMD quintile with single morbidity\n booster vaccination and no previous COVD")
ggplot2::ggsave( filename = here::here("output","dailyVTERisk.pdf"),device = "pdf",width = 8, height = 8, units = 'in',dpi = 'retina')




