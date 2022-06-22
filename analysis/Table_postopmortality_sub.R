library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

load(file = here::here("output","cohort_long.RData"))
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

post.op.died.model.sub <- 
  list(survival::coxph(survival::Surv(start,end,died) ~ Colectomy + Cholecystectomy + HipReplacement + KneeReplacement + 
                         postcovid*wave + age.cat + sex + bmi.cat + imd5 + vaccination.status.factor + region + Current.Cancer +
                         Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                       data = dt.tv[start >=0 & sub.op == T & any.op == T], model = T))
data.table::fwrite(broom::tidy(post.op.died.model.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopdiedmodelsub.csv"))


new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(30,8*length(procedures)),
                                                'event' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
                                                'postcovid' = rep(c(rep(F,4),rep(T,4)),
                                                                  times = length(procedures)),
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

n.type.events <- 1
cuminc.adjusted.mortality.sub <-   matrix(cuminc.cox(n.type.events = n.type.events,
                                                 dt = 'dt.tv[start >=0  & sub.op == T & any.op == T]',
                                                 model = 'post.op.died.model.sub',
                                                 newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)


colnames(cuminc.adjusted.mortality.sub) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.mortality.sub) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.died.model.sub,cuminc.adjusted.mortality.sub, file = here::here("output","postopmortality_sub.RData"))