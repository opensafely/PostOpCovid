library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))
set.seed(23097859)
###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures.sub <- c('Colectomy','Cholecystectomy',
                    'HipReplacement','KneeReplacement')
covariates <- c(procedures.sub,'sex','age.cat','bmi.cat','imd5','wave','LOS.bin',
                'vaccination.status.factor','Current.Cancer','Emergency','Charl12','recentCOVID','previousCOVID','region')

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','postop.covid.cohort','end.fu'))]

dt.tv[,(drop.vars) := NULL]


dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))

dt.tv[,(procedures.sub) := lapply(.SD,function(x) x==1), .SDcols = (procedures.sub)]
dt.tv <- dt.tv[(postop.covid.cohort) & sub.op == T]
gc()

n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) & sub.op == T ,event]))[-1]

library(doParallel)
ncores <- parallel::detectCores(logical = F)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


samples.sub <-  foreach::foreach(
  i = 1:1000,
  .combine = cbind,
  .multicombine = T,
  .inorder = F,
  .verbose = F,
  .packages = c('data.table', 'survival','foreach'),
  .export = c(
    'n.type.events',
    'dt.tv',
    'procedures.sub',
    'covariates'
  )
) %dopar% {
  dt.tv <-
    dt.tv[patient_id %in% sample(unique(patient_id), replace = T) &
            (postop.covid.cohort)]
  post.op.covid.model.sub <-
    lapply(n.type.events, function(i)
      survival::coxph(
        survival::Surv(start, end, event == i) ~ Colectomy + Cholecystectomy  + KneeReplacement +
          age.cat + sex + bmi.cat + imd5 +
          vaccination.status.factor + region + Current.Cancer +
          Emergency + LOS.bin + wave + Charl12 + recentCOVID + previousCOVID,
        id = patient_id,
        data = dt.tv[(postop.covid.cohort)],
        model = T
      ))
  data.table::as.data.table(foreach::foreach(
    predi = 1:length(covariates),
    .combine = 'rbind',
    .inorder = T
  ) %do% {

    newdata.rows <-
      length(unique(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])]))
    
    
    newdata.pred <-
      data.table::data.table(
        'start' = rep(0, newdata.rows),
        'end' = rep(30, newdata.rows),
        'event' = rep(F, newdata.rows),
        'patient_id' = 1:newdata.rows,
        'Colectomy' = c(rep(T,newdata.rows)),
        'Cholecystectomy'=c(rep(F,newdata.rows)),
        'HipReplacement'=c(rep(F,newdata.rows)),
        'KneeReplacement'=c(rep(F,newdata.rows)),
        'age.cat' = rep('(50,70]', newdata.rows),
        'sex' = rep('F', newdata.rows),
        'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2], newdata.rows),
        'imd5' = rep(levels(dt.tv$imd5)[3], newdata.rows),
        'wave' = rep(paste0('Wave_', 3), times = newdata.rows),
        'vaccination.status.factor' = rep('3', newdata.rows),
        'region' = rep("East Midlands", newdata.rows),
        'Current.Cancer' = rep(T, newdata.rows),
        'LOS.bin' = rep(F, newdata.rows),
        'Emergency' =  rep(F, newdata.rows),
        'Charl12' =  rep('Single', newdata.rows),
        'recentCOVID' = rep(F, newdata.rows),
        'previousCOVID' = rep(F, newdata.rows)
      )
    if (predi <= length(procedures.sub)) {
      newdata.pred[, (procedures.sub) := F]
      newdata.pred[, (procedures.sub[predi]) := c(F, T)]
    } else {
      if (is.factor(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])])) {
        newdata.pred[, (covariates[predi]) :=  as.character(sort(unique(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])], na.rm = T)))]
      } else if (is.logical(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])])) {
        newdata.pred[, (covariates[predi]) :=  as.logical(sort(unique(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])], na.rm = T)))]
      } else if (is.numeric(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])])) {
        newdata.pred[, (covariates[predi]) :=  is.numeric(sort(unique(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])], na.rm = T)))]
      } else {
        newdata.pred[, (covariates[predi]) := sort(unique(dt.tv[!is.na(get(covariates[predi])), get(covariates[predi])], na.rm = T))]
      }
    }
    
    cuminc.cox(
      n.type.events = n.type.events,
      dt = 'dt.tv',
      model = 'post.op.covid.model.sub',
      newdata = 'newdata.pred',
      day = 30
    )
  
  })
}

t.samples <- t(apply(samples.sub,1,quantile,c(0.25,0.5,0.75)))
boot.IQR.sub <-apply(t.samples,1,function(x) paste0(x[2],' (',x[1],',',x[3],')'))


save(samples.sub, file = here::here("output","postopcovid_adjusted_bootstrappedIQR_sub.RData"))

data.table::fwrite(boot.IQR.sub, file = here::here("output","postopcovid_adjusted_bootstrappedIQR_sub.csv"))
