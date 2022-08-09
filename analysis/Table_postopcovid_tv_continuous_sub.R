library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

#dt.tv.splits <- data.table::setDT(arrow::read_feather(here::here("output","cohort_postdisch_week_splits.feather")))
#dt.tv.splits <- dt.tv.splits[,.(patient_id,Colectomy,Cholecystectomy,HipReplacement,KneeReplacement,age.cat,sex,bmi.cat,imd5,wave,vaccination.status.factor,region,Current.Cancer,Emergency,
                                #week.post.disch,week.post.op,LOS.bin,Charl12,recentCOVID,previousCOVID,postcovid,start.readmit,end.readmit,
                               # tstart,tstop,end.fu,start,end,event,postop.covid.cohort,los.end,event.VTE,postcovid.VTE.cohort,study.start)]
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')
days.limit <- 35

# Start = 0 = day of operation
#dt.tv.splits[, `:=`(start = tstart - study.start,
                    #end = tstop - study.start)]
#dt.tv.splits <- dt.tv.splits[start >= 0,] # Need to start follow up on day after operation as can't identify order when events on same day

dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
               (is.finite(Cholecystectomy) & Cholecystectomy == T) |
               (is.finite(HipReplacement)  & HipReplacement == T) | 
               (is.finite(KneeReplacement) & KneeReplacement == T)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))

#dt.tv.splits[event == 3, event := 2]
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) & sub.op == T,event]))[-1]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

data.table::setkey(dt.tv, patient_id, end.fu, start)
post.op.covid.model.split.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event==i) ~  Colectomy + Cholecystectomy +  
                                                      KneeReplacement + age.cat + 
  sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency +  Charl12 + recentCOVID + previousCOVID,
   id = patient_id,
  data = dt.tv[(postop.covid.cohort) & start <=end & sub.op == T], model = T))

names(post.op.covid.model.split.sub) <- c('COVID-19','Non COVID-19 emergency readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.covid.model.split.sub,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, conf_level = .95, exponentiate = TRUE, output = here::here("output","postopcovidmodeldailysub.html"))

newdata.rows <- 1 #

newdata.pred <- data.table::data.table('start'  = 0,
                                       'end' = days.limit,
                                       'event' = rep(F,newdata.rows),
                                   #   'week.post.op' = paste(rep(1:5,each = 7)),
                                      'patient_id' = 1:newdata.rows,
                                      'Colectomy' =  rep(T,newdata.rows),
                                      'Cholecystectomy'= rep(F,newdata.rows),
                                      'HipReplacement'= rep(F,newdata.rows),
                                      'KneeReplacement'= rep(F,newdata.rows),
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
# 
# base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.covid.model.split.sub[[i]],centered = F))
# 
# base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
#                                                                             base.haz = base.haz[[i]][,1] - 
#                                                                               c(0,head(base.haz[[i]][,1],-1)))})
# 
# lp <- lapply(n.type.events, function(i) {  data.table::data.table('time' = 0:34,
#                                                                     'risk' = exp(predict(object = post.op.covid.model.split.sub[[i]],
#                                                                                          type = 'lp', 
#                                                                                          newdata = newdata.pred))) })
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[[j]])), j, 0)
# 
# daily.post.op.risk.sub <- 
#   unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
#     lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] 
#   })),2,cumsum))*
#     lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,2]),.SDcols = c(2:3)] ,2,cumsum), digits = 3))

daily.post.op.risk.sub <- cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postop.covid.cohort) & start <=end & sub.op == T]', model = 'post.op.covid.model.split.sub', newdata = 'newdata.pred', day = 1:days.limit)

times.comb <-sort(unique(dt.tv[(postop.covid.cohort) & start <=end & sub.op == T][ event %in% n.type.events & end <= days.limit ,end]))

daily.post.op.risk.sub[!is.finite(daily.post.op.risk.sub)] <- 0


daily.post.op.risk.sub  <-  data.table::data.table(`risk` = daily.post.op.risk.sub - c(0,daily.post.op.risk.sub[-length(daily.post.op.risk.sub)]),
                                                `Days post op` = times.comb)[times.comb <=days.limit,]

data.table::fwrite(daily.post.op.risk.sub, file = here::here("output","daily_postopcovid_tv_sub.csv"))

ggplot2::ggplot(daily.post.op.risk.sub[`Days post op`<= days.limit]) +
  ggplot2::geom_line(ggplot2::aes(x = `Days post op`, y = risk)) + 
  ggplot2::geom_smooth(ggplot2::aes(x = `Days post op`, y = risk)) +
  ggplot2::ylab("Daily risk (%)") +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle("Daily risk of COVID from day of operation",
                   subtitle = "Predicted for  elective colectomy cancer operation in 50-70 year old female\n BMI 20-25, 3rd IMD quintile with single morbidity\n booster vaccination and no previous COVD")
ggplot2::ggsave( filename = here::here("output","dailyCovidRiskSub.pdf"),device = "pdf",width = 8, height = 8, units = 'in',dpi = 'retina')

##############
# Not enough deaths to treat separately from emergency readmissions
dt.tv.splits[event.VTE == 3, event.VTE := 2]
n.type.events <- sort(unique(dt.tv[(postcovid.VTE.cohort) & sub.op == T,event.VTE]))[-1]


post.op.VTE.model.split.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.VTE==i) ~  Colectomy + postcovid +    
                                                      Cholecystectomy + KneeReplacement + age.cat + 
                                            sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + 
                                            Current.Cancer + Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID, 
                                          id = patient_id,
                                          data = dt.tv[(postcovid.VTE.cohort) & start <=end], model = T))

names(post.op.VTE.model.split.sub) <- c('Post discharge VTE','Non COVID-19 emergency readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.VTE.model.split.sub,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodeldailysub.html"))


newdata.rows = 1

newdata.pred <- data.table::data.table('start.readmit' = rep(0, times = 2),
                                       'end.readmit' = rep(days.limit,times = 2),
                                       'event.VTE' = rep(F,newdata.rows*2),
                               #        'week.post.disch' = paste(rep(ceiling(1:newdata.rows / 7), times = 2)),
                                       'patient_id' = rep(1:2,each = newdata.rows),
                                       'Colectomy' =  rep(T,newdata.rows*2),
                                       'Cholecystectomy'= rep(F,newdata.rows*2),
                                       'HipReplacement'= rep(F,newdata.rows*2),
                                       'KneeReplacement'= rep(F,newdata.rows*2),
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
# 
# base.haz <- lapply(n.type.events, function(i) survival::basehaz(post.op.VTE.model.split.sub[[i]],centered = F))
# 
# base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
#                                                                             base.haz = base.haz[[i]][,1] - 
#                                                                               c(0,head(base.haz[[i]][,1],-1)))})
# 
# lp <- lapply(n.type.events, function(i) {  data.table::dcast(
#   data.table::data.table('patient_id' = rep(1:2, each = newdata.rows),
#   'time' = rep(0:(newdata.rows - 1),times = 2),
#   'risk' = exp(predict(object = post.op.VTE.model.split.sub[[i]],
#                         type = 'lp', 
#                         newdata = newdata.pred))),time ~patient_id, value.var = 'risk')})
# 
# base.haz.merge <- Reduce(x =base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T, sort = T))
# 
# for (j in 1:ncol(base.haz.merge)) set(base.haz.merge, which(!is.finite(base.haz.merge[,(j)])), j, 0)
# 
# daily.post.op.VTE.risk.sub <- 
#   unlist(round(100*apply(exp(apply(-Reduce('+',lapply(n.type.events, function(i) {
#     lp[[i]][base.haz.merge[order(time),.SD,.SDcols = c(1,i+1)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)] 
#   })),2,cumsum))*
#     lp[[1]][base.haz.merge[order(time),.SD,.SDcols = c(1,2)],,roll =Inf,on = 'time', rollends = c(T,T)][time >= 0][order(time),.(.SD[,1]*.SD[,3], .SD[,2]*.SD[,3]),.SDcols = c(2:4)], 2, cumsum ), digits = 3))

daily.post.op.VTE.risk.sub <- cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & start <=end & sub.op == T]', model = 'post.op.VTE.model.split.sub', newdata = 'newdata.pred', day = 1:days.limit)

times.comb <-sort(unique(dt.tv[(postcovid.VTE.cohort) & start <=end & sub.op == T][ event.VTE %in% n.type.events & end <=days.limit,end]))

daily.post.op.VTE.risk.sub[!is.finite(daily.post.op.VTE.risk.sub)] <- 0


daily.post.op.VTE.risk.sub <- as.vector(daily.post.op.VTE.risk.sub - rbind(c(0,0),daily.post.op.VTE.risk.sub[-nrow(daily.post.op.VTE.risk.sub),]))
daily.post.op.VTE.risk.sub  <-  cbind(data.table::data.table("COVID"= rep(c(F,T), each = length(times.comb)),
"Risk" = daily.post.op.VTE.risk.sub,
"Days post discharge" = rep(times.comb, times = 2)))

##################################
save(daily.post.op.risk.sub,daily.post.op.VTE.risk.sub, file = here::here("output","postopcovid_tv_daily_sub.RData"))
data.table::fwrite(daily.post.op.VTE.risk.sub, file = here::here("output","daily_postopcovid_VTE_tv_sub.csv"))

######################
ggplot2::ggplot(daily.post.op.VTE.risk.sub[`Days post discharge`<= newdata.rows]) +
  ggplot2::geom_line(ggplot2::aes(x = `Days post discharge`, y = Risk, colour = COVID)) + 
  ggplot2::geom_smooth(ggplot2::aes(x = `Days post discharge`, y = Risk, colour = COVID)) +
  ggplot2::ylab("Daily risk (%)") +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle("Daily risk of VTE from day of discharge", 
  subtitle = "Predicted for elective colectomy cancer operation in 50-70 year old female\n BMI 20-25, 3rd IMD quintile with single morbidity\n booster vaccination and no previous COVD")
ggplot2::ggsave( filename = here::here("output","dailyVTERiskSub.pdf"),device = "pdf",width = 8, height = 8, units = 'in',dpi = 'retina')


