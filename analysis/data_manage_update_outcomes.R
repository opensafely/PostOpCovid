
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))

index_date <- data.table::as.IDate("2020-02-01")
last_date <- data.table::as.IDate("2022-03-01")

procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

dt.update.outcomes <- Reduce(function(...) {
  merge(..., by = c('patient_id'), all = T)
}, lapply(procs, function(x) data.table::fread(here::here('output',
                                                          paste0('input_',x,'.csv')))[region !='',][, c('date_admitted','date_discharged',
                                                                                                        'dob','bmi_date_measured',
                                                                                                        'date_death_ons',
                                                                                                        'date_death_cpns',
                                                                                                        'age','region','sex',
                                                                                                        'bmi','imd','died') := NULL]))


data.table::setkey(dt.update.outcomes,patient_id)
arrow::write_feather(dt.update.outcomes, sink = here::here("output","update_outcomes.feather"))
dt <- data.table('patient_id' = dt.update.outcomes[,patient_id])
rm(dt.update.outcomes)

dt.major <- Reduce(function(...) {
  merge(..., by = c('patient_id'), all = T)
},lapply(procs,
         function(x) tryCatch(data.table::fread(here::here('output',
                                                  paste0('input_',x,'_majorminor.csv')))[,
                                                  c('date_admitted','date_discharged') := NULL],
                                                  error = function(e) EVAL(paste0("dt[,.('patient_id' = patient_id,'",x,"_Major_HES_binary_flag'  = NA)]")))
                                                  ))


data.table::setkey(dt.major,patient_id)

arrow::write_feather(dt.major, sink = here::here("output","update_major.feather"))
rm(dt.major)