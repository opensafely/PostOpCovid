library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(2)

source(here::here("analysis","Utils.R"))
set.seed(23097859)
###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')


covariates <- c(procedures,'sex','age.cat','bmi.cat','imd5','wave',
                'vaccination.status.factor','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID','region')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','postop.covid.cohort','end.fu'))]

dt.tv[,(drop.vars) := NULL]
dt.tv <- dt.tv[(postop.covid.cohort)]
gc()
n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) ,event]))[-1]

library(doParallel)
ncores <- parallel::detectCores(logical = F)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


samples <-  foreach::foreach(
  i = 1:1000,
  .combine = cbind,
  .multicombine = T,
  .inorder = F,
  .verbose = F,
  .packages = c('data.table', 'survival','foreach'),
  .export = c(
    'n.type.events',
    'dt.tv',
    'procedures',
    'covariates'
  )
) %dopar% {
  dt.tv <-
    dt.tv[patient_id %in% sample(unique(patient_id), replace = T) &
            (postop.covid.cohort)]
  post.op.covid.model <-
    lapply(n.type.events, function(i)
      survival::coxph(
        survival::Surv(start, end, event == i) ~ Abdominal + Cardiac +
          Obstetrics + Thoracic + Vascular +
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
        'Abdominal' = rep(T, newdata.rows),
        'Cardiac' = rep(F, newdata.rows),
        'Obstetrics' =
          rep(F, newdata.rows),
        'Orthopaedic' =
          rep(F, newdata.rows),
        'Thoracic' =
          rep(F, newdata.rows),
        'Vascular' =
          rep(F, newdata.rows),
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
    if (predi <= length(procedures)) {
      newdata.pred[, (procedures) := F]
      newdata.pred[, (procedures[predi]) := c(F, T)]
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
    
    result <- try(cuminc.cox(
      n.type.events = n.type.events,
      dt = 'dt.tv',
      model = 'post.op.covid.model',
      newdata = 'newdata.pred',
      day = 30
    ), silent = T)

    if(class(result)[1] == "try-error") {
     return(rep(NA,newdata.rows))
    } else {
      return(result)
    }
  
  })
}

t.samples <- t(apply(samples,1,quantile,c(0.25,0.5,0.75), na.rm = T))
boot.IQR <-apply(t.samples,1,function(x) paste0(x[2],' (',x[1],',',x[3],')'))


save(samples, file = here::here("output","postopcovid_adjusted_bootstrappedIQR.RData"))

data.table::fwrite(boot.IQR, file = here::here("output","postopcovid_adjusted_bootstrappedIQR.csv"))
