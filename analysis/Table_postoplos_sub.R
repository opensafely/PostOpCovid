library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)
library(flexsurv)
source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

covariates <- c(procedures,'age.cat','sex','bmi.cat','imd5','wave',
                'vaccination.status.factor','region','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))


dt.tv[,event.los := (died == T)*2]
dt.tv[event.los !=2,event.los := discharged == T]
dt.tv[, postop.los.cohort := start>=0 & tstop <= los.end & end <= 90 & any.op ==T]

n.type.events <- sort(unique(dt.tv[(postop.los.cohort) ,event.los]))[-1]

post.op.LOS.model.sub <-  flexsurv::flexsurvreg(survival::Surv(start,end, event.los == 1) ~ Colectomy + Cholecystectomy +  KneeReplacement + 
                                              postcovid + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + 
                                              Emergency + Charl12 + recentCOVID + previousCOVID, 
                                           data = dt.tv[(postop.los.cohort) & sub.op == T],
                                           dist = 'weibull')

#  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end, event.los == i) ~ Abdominal + Cardiac + Obstetrics + Orthopaedic + Thoracic + Vascular + postcovid + age.cat + sex + bmi.cat + imd5 + wave + vaccination.status.factor + region + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,data = dt.tv[(postop.los.cohort)], model = T))


#data.table::fwrite(broom::tidy(post.op.LOS.model, exponentiate= T, conf.int = T), file = here::here("output","post.op.LOS.model.csv"))
write.table(tidy.flexsurvreg(post.op.LOS.model.sub), here::here("output","postopLOSmodel_sub.csv"), sep = ",",quote = F, row.names = F)

new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(30,8*length(procedures)),
                                                'event' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(1,8),rep(0,24)),
                                                'Cholecystectomy'=c(rep(0,8),rep(1,8),rep(0,16)),
                                                'HipReplacement'=c(rep(0,16),rep(1,8),rep(0,8)),
                                                'KneeReplacement'=c(rep(0,24),rep(1,8)),
                                                'postcovid' = as.numeric(rep(c(rep(F,4),rep(T,4)), times = length(procedures))),
                                                'age.cat' = rep('(50,70]',8*length(procedures)),
                                                'sex' = rep('F',8*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                'region' = rep("East Midlands",8*length(procedures)),
                                                'Current.Cancer' = rep(T,8*length(procedures)),
                                                'Emergency' =   rep(F,8*length(procedures)),
                                                'Charl12' =  rep('Single',8*length(procedures)),
                                                'recentCOVID' = rep(F,8*length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))



mean.adjusted.los.sub <- summary(object = post.op.LOS.model.sub, newdata = new.data.postop.covid, type = "mean", ci = T, tidy = T)[,1:3]
 # matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv', model = 'post.op.LOS.model', newdata = 'new.data.postop.covid', day = 14), byrow = T, ncol = 4)

mean.adjusted.los.sub <- matrix(apply(mean.adjusted.los.sub,1,function(x) paste0(round(x[1], digits = 1), 
                                                                  " days, (",
                                                                  round(x[2], digits = 1),
                                                                  ',',
                                                                  round(x[3], digits = 1),')')), byrow = T, ncol = 4)

colnames(mean.adjusted.los.sub) <- paste0('Wave_',1:4)
rownames(mean.adjusted.los.sub) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))
save(mean.adjusted.los.sub,post.op.LOS.model.sub, file = here::here("output","postoplos_sub.RData"))
data.table::fwrite(mean.adjusted.los.sub, file = here::here("output","postoplos_sub.csv"))
