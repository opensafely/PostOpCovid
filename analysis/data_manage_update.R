
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))

procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

dt.update <- Reduce(function(...) {
  merge(..., by = c('patient_id'), all = T)
  },lapply(procs, 
           function(x) data.table::fread(here::here('output',
                                                    paste0('input_',x,'.csv')))[,(paste0(x,'_date_admitted')) := date_admitted][,
                                                                                 (paste0(x,'_date_discharged')) := date_discharged][
                                                                   region !='',][,                                                                
                                                                    c('date_admitted','date_discharged',
                                                                    'dob','bmi_date_measured',
                                                                    'date_death_ons',
                                                                    'date_death_cpns',
                                                                    'age','region','sex',
                                                                    'bmi','imd','died') := NULL]))

data.table::setkey(dt.update,patient_id)

summary(dt.update)

## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.update)[grepl(pattern = paste(procedures,collapse = '|'),x = names(dt.update))]
#non.op.vars <- names(dt.update)[!grepl(pattern = paste(procedures,collapse = '|'),x = names(dt.update))]

dt.update <- data.table::melt(dt.update, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures,'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),collapse = "|"),
                                                                                        replacement = "",
                                                                                        x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures,unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),
                                                                                                              collapse = "|"),
                                                                                             replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]
#dt.update <- dt.update[,non.op.vars, with = F][aggregate_operations, on = "patient_id"]
#rm(aggregate_operations)




## Baseline fixed factors ----
fixed <- c('patient_id')

## Start dates of exposures associated with each admission: needs to start from beginning of row time period ----

proc.time.stubs.start <- c('_date_admitted',
                           '_recent_date',
                           '_previous_date')
proc.time.cols.start <- paste(rep(procedures,each = length(proc.time.stubs.start)),proc.time.stubs.start, sep ="")

## Outcome dates so need to be flagged at end of row time period ----
proc.time.stubs.end <- c('_date_discharged',
                         '_emergency_readmit_date_admitted',
                         '_date',
                         '_VTE_HES_date_admitted',
                         '_VTE_GP_date',
                         '_anticoagulation_prescriptions_date')
proc.time.cols.end <- paste(rep(procedures,each = length(proc.time.stubs.end)),proc.time.stubs.end, sep ="")

# Set dates as numeric to avoid any problems with comparisons (was an issue for tmerge, probably not important with rolling dates)
dt.update[,(c(proc.time.cols.start,proc.time.cols.end)) := lapply(.SD,as.numeric), .SDcols = c(proc.time.cols.start,proc.time.cols.end)]

## Final date for per patient taking into account final patient in GP database and currrent study end date ----

max.grp.col_(dt = 'dt.update', max.var.name = 'max.date', aggregate.cols = c(proc.time.cols.end,proc.time.cols.start), id.vars = c('patient_id')) 
min.grp.col_(dt = 'dt.update', min.var.name = 'min.date', aggregate.cols = c(proc.time.cols.end,proc.time.cols.start), id.vars = c('patient_id')) 

# maximum date of follow up in data 
max.date.fu <- max(as.numeric(data.table::as.IDate(dt.update$max.date)), na.rm = T)

dt.update[, max.date := do.call(pmax, c(.SD, na.rm = T)), .SDcols = proc.time.cols.end]
dt.update[!is.finite(max.date), max.date := as.numeric(data.table::as.IDate('2022-03-01') + 90)]
dt.update[,tstart := do.call(pmin, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_admitted")]
dt.update[op.number == 1 ,tstart:= min.date]

# maximum date of follow up in data 
max.date.fu <- max(as.numeric(dt.update$max.date), na.rm = T)


## All times for long cohort table ----
dt.dates.wide <- dt.update[,.SD, .SDcols = c('patient_id',
                                      proc.time.cols.start, 
                                      proc.time.cols.end)]



# Roll dates into main data ----

# melt dates to long
dt.dates.long <- unique(data.table::melt(dt.dates.wide, 
                                         id.vars = 'patient_id',na.rm = T, 
                                         value.name = 'tstart')[is.finite(tstart),c(1,3)])
data.table::setkey(data.table::setDT(dt.update),patient_id, tstart)
data.table::setkey(dt.dates.long, patient_id, tstart)

# Roll dates into data with locf for all variables, (and nocb for any earlier times)
dt.tv.update <- dt.update[dt.dates.long,,
                          rollends = c(T,T), 
                          roll = Inf, 
                          on = c('patient_id','tstart'),
                          mult = 'all'][,
                          c('max.date','min.date',"op.number") := NULL]

arrow::write_feather(dt.tv.update, sink = here::here("output","cohort_long_update.feather"))

