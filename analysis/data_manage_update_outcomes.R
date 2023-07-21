
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))

index_date <- data.table::as.IDate("2020-02-01")
last_date <- data.table::as.IDate("2022-03-01")

procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

dt.update.outcomes <- 
  Reduce(function(...) {
    merge(..., by = c('patient_id'), all = T)}, Filter(Negate(is.null),lapply(procs, function(p) {
      temp <- read.csv(file = here::here('output',
                                         paste0('input_',p,'.csv')),header = T,sep = ',',fill = T)
      data.table::setDT(temp)
      temp[,c('date_admitted','date_discharged','dob','bmi_date_measured',
              'date_death_ons',
              'date_death_cpns',
              'age','sex',
              'bmi','imd','died'):= NULL]
      temp <- temp[region !='',]
      temp[,region := NULL]
      date.cols <- names(temp)[grepl('date',names(temp),ignore.case = T)]
      temp[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]
      return(temp)
      rm(temp)  
    })) 
   )

data.table::setkey(dt.update.outcomes,patient_id)
arrow::write_feather(dt.update.outcomes, sink = here::here("output","update_outcomes.feather"))
dt <- data.table('patient_id' = dt.update.outcomes[,patient_id])
rm(dt.update.outcomes)


dt.major <- 
  Reduce(function(...) {
    merge(..., by = c('patient_id'), all = T)}, Filter(Negate(is.null),lapply(procs, function(p) {
      temp <- tryCatch(read.csv(file = here::here('output',
                       paste0('input_',p,'_majorminor.csv')),header = T,sep = ',',fill = T),
                       error = function(e) EVAL(paste0("dt[,.('patient_id' = patient_id,'",p,"_Major_HES_binary_flag'  = NA,date_admitted = NA, date_discharged = NA)]")))
      data.table::setDT(temp)
      temp[,c('date_admitted','date_discharged'):= NULL]
      return(temp)
      rm(temp)  
    })) 
  )

data.table::setkey(dt.major,patient_id)

arrow::write_feather(dt.major, sink = here::here("output","update_major.feather"))
rm(dt.major)

procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

#### Abdomen
dt.Abdo <- read.csv(file = here::here('output',
                                        paste0('input_Abdominal.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Abdo)
dt.Abdo <- dt.Abdo[region !='',]
dt.Abdo[,region := NULL]

date.cols <- names(dt.Abdo)[grepl('date',names(dt.Abdo),ignore.case = T)]
dt.Abdo[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]

## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[1],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[1],paste0("_",1:5),paste0),'_date_discharged'))
for (i in 1:length(op.admit.vars)) {
  dt.Abdo[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
          (names(dt.Abdo)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                      pattern = '_date_admitted',
                                                      replacement = ""),'*'),
                                x = names(dt.Abdo))]) := NA]
}

## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Abdo)[grepl(pattern = paste(procedures[1],collapse = '|'),x = names(dt.Abdo))]
non.op.vars <- names(dt.Abdo)[!grepl(pattern = paste(procedures[1],collapse = '|'),x = names(dt.Abdo))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Abdo, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[1],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[1],unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]



vars_char <- names(dt.Abdo[,non.op.vars,with = F])[sapply(dt.Abdo[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Abdo[,non.op.vars,with = F])[sapply(dt.Abdo[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Abdo[,non.op.vars,with = F])[sapply(dt.Abdo[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Abdo[,non.op.vars,with = F])[sapply(dt.Abdo[,non.op.vars,with = F], is.double)]


dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Abdo, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Abdo, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Abdo, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Abdo, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Abdo <- dt.baseline[aggregate_operations, on = "patient_id"]

#### Cardiac
dt.Cardiac <- read.csv(here::here('output',
                                           paste0('input_Cardiac.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Cardiac)
dt.Cardiac <- dt.Cardiac[region !='',]
dt.Cardiac[,region := NULL]

date.cols <- names(dt.Cardiac)[grepl('date',names(dt.Cardiac),ignore.case = T)]
dt.Cardiac[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]


## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[2],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[2],paste0("_",1:5),paste0),'_date_discharged'))

for (i in 1:length(op.admit.vars)) {
  dt.Cardiac[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
             (names(dt.Cardiac)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                            pattern = '_date_admitted',
                                                            replacement = ""),'*'),
                                      x = names(dt.Cardiac))]) := NA]
}
## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Cardiac)[grepl(pattern = paste(procedures[2],collapse = '|'),x = names(dt.Cardiac))]
non.op.vars <- names(dt.Cardiac)[!grepl(pattern = paste(procedures[2],collapse = '|'),x = names(dt.Cardiac))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Cardiac, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[2],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures[2],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[2],unique(gsub(pattern = paste0(as.vector(outer(procedures[2],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]


vars_char <- names(dt.Cardiac[,non.op.vars,with = F])[sapply(dt.Cardiac[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Cardiac[,non.op.vars,with = F])[sapply(dt.Cardiac[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Cardiac[,non.op.vars,with = F])[sapply(dt.Cardiac[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Cardiac[,non.op.vars,with = F])[sapply(dt.Cardiac[,non.op.vars,with = F], is.double)]

dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Cardiac, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Cardiac, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Cardiac, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Cardiac, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Cardiac <- dt.baseline[aggregate_operations, on = "patient_id"]

#### Obstetrics
dt.Obstetrics<- read.csv(here::here('output',
                                             paste0('input_Obstetrics.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Obstetrics)
dt.Obstetrics <- dt.Obstetrics[region !='',]
dt.Obstetrics[,region := NULL]

date.cols <- names(dt.Obstetrics)[grepl('date',names(dt.Obstetrics),ignore.case = T)]
dt.Obstetrics[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]

## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[3],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[3],paste0("_",1:5),paste0),'_date_discharged'))

for (i in 1:length(op.admit.vars)) {
  dt.Obstetrics[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
                (names(dt.Obstetrics)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                                  pattern = '_date_admitted',
                                                                  replacement = ""),'*'),
                                            x = names(dt.Obstetrics))]) := NA]
}
## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Obstetrics)[grepl(pattern = paste(procedures[3],collapse = '|'),x = names(dt.Obstetrics))]
non.op.vars <- names(dt.Obstetrics)[!grepl(pattern = paste(procedures[3],collapse = '|'),x = names(dt.Obstetrics))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Obstetrics, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[3],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures[3],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[3],unique(gsub(pattern = paste0(as.vector(outer(procedures[3],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]


vars_char <- names(dt.Obstetrics[,non.op.vars,with = F])[sapply(dt.Obstetrics[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Obstetrics[,non.op.vars,with = F])[sapply(dt.Obstetrics[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Obstetrics[,non.op.vars,with = F])[sapply(dt.Obstetrics[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Obstetrics[,non.op.vars,with = F])[sapply(dt.Obstetrics[,non.op.vars,with = F], is.double)]

dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Obstetrics, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Obstetrics, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Obstetrics, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Obstetrics, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Obstetrics <- dt.baseline[aggregate_operations, on = "patient_id"]

#### Orthopaedic
dt.Orthopaedic <- read.csv(here::here('output',
                                               paste0('input_Orthopaedic.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Orthopaedic)
dt.Orthopaedic <- dt.Orthopaedic[region !='',]
dt.Orthopaedic[,region := NULL]

date.cols <- names(dt.Orthopaedic)[grepl('date',names(dt.Orthopaedic),ignore.case = T)]
dt.Orthopaedic[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]


## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[4],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[4],paste0("_",1:5),paste0),'_date_discharged'))

for (i in 1:length(op.admit.vars)) {
  dt.Orthopaedic[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
                 (names(dt.Orthopaedic)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                                    pattern = '_date_admitted',
                                                                    replacement = ""),'*'),
                                              x = names(dt.Orthopaedic))]) := NA]
}
## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Orthopaedic)[grepl(pattern = paste(procedures[4],collapse = '|'),x = names(dt.Orthopaedic))]
non.op.vars <- names(dt.Orthopaedic)[!grepl(pattern = paste(procedures[4],collapse = '|'),x = names(dt.Orthopaedic))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Orthopaedic, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[4],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures[4],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[4],unique(gsub(pattern = paste0(as.vector(outer(procedures[4],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]


vars_char <- names(dt.Orthopaedic[,non.op.vars,with = F])[sapply(dt.Orthopaedic[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Orthopaedic[,non.op.vars,with = F])[sapply(dt.Orthopaedic[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Orthopaedic[,non.op.vars,with = F])[sapply(dt.Orthopaedic[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Orthopaedic[,non.op.vars,with = F])[sapply(dt.Orthopaedic[,non.op.vars,with = F], is.double)]

dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Orthopaedic, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Orthopaedic, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Orthopaedic, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Orthopaedic, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Orthopaedic <- dt.baseline[aggregate_operations, on = "patient_id"]

#### Thoracic
dt.Thoracic <- read.csv(here::here('output',
                                            paste0('input_Thoracic.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Thoracic)
dt.Thoracic <- dt.Thoracic[region !='',]
dt.Thoracic[,region := NULL]

date.cols <- names(dt.Thoracic)[grepl('date',names(dt.Thoracic),ignore.case = T)]
dt.Thoracic[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]


## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[5],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[5],paste0("_",1:5),paste0),'_date_discharged'))

for (i in 1:length(op.admit.vars)) {
  dt.Thoracic[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
              (names(dt.Thoracic)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                              pattern = '_date_admitted',
                                                              replacement = ""),'*'),
                                        x = names(dt.Thoracic))]) := NA]
}
## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Thoracic)[grepl(pattern = paste(procedures[5],collapse = '|'),x = names(dt.Thoracic))]
non.op.vars <- names(dt.Thoracic)[!grepl(pattern = paste(procedures[5],collapse = '|'),x = names(dt.Thoracic))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Thoracic, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[5],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures[5],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[5],unique(gsub(pattern = paste0(as.vector(outer(procedures[5],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]


vars_char <- names(dt.Thoracic[,non.op.vars,with = F])[sapply(dt.Thoracic[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Thoracic[,non.op.vars,with = F])[sapply(dt.Thoracic[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Thoracic[,non.op.vars,with = F])[sapply(dt.Thoracic[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Thoracic[,non.op.vars,with = F])[sapply(dt.Thoracic[,non.op.vars,with = F], is.double)]

dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Thoracic, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Thoracic, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Thoracic, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Thoracic, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Thoracic <- dt.baseline[aggregate_operations, on = "patient_id"]

#### Vascular
dt.Vascular <- read.csv(here::here('output',
                                            paste0('input_Vascular.csv')),header = T,sep = ',',fill = T)
data.table::setDT(dt.Vascular)
dt.Vascular <- dt.Vascular[region !='',]
dt.Vascular[,region := NULL]

date.cols <- names(dt.Vascular)[grepl('date',names(dt.Vascular),ignore.case = T)]
dt.Vascular[,(date.cols) := lapply(.SD,function(x) as.IDate(x,'%Y-%m-%d')), .SDcols = date.cols]


## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures[6],paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures[6],paste0("_",1:5),paste0),'_date_discharged'))

for (i in 1:length(op.admit.vars)) {
  dt.Vascular[get(op.admit.vars[i])>get(op.discharge.vars[i]) | data.table::as.IDate(get(op.admit.vars[i])) > data.table::as.IDate('2022-10-01'),
              (names(dt.Vascular)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                                              pattern = '_date_admitted',
                                                              replacement = ""),'*'),
                                        x = names(dt.Vascular))]) := NA]
}

## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt.Vascular)[grepl(pattern = paste(procedures[6],collapse = '|'),x = names(dt.Vascular))]
non.op.vars <- names(dt.Vascular)[!grepl(pattern = paste(procedures[6],collapse = '|'),x = names(dt.Vascular))]
non.op.vars <- non.op.vars[!(non.op.vars %in% 'patient_id')]

aggregate_operations <- data.table::melt(dt.Vascular, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures[6],'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures[6],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures[6],unique(gsub(pattern = paste0(as.vector(outer(procedures[6],paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]


vars_char <- names(dt.Vascular[,non.op.vars,with = F])[sapply(dt.Vascular[,non.op.vars,with = F], is.character)]
vars_date <- names(dt.Vascular[,non.op.vars,with = F])[sapply(dt.Vascular[,non.op.vars,with = F], function(x) inherits(x, 'Date') )]
vars_int <- names(dt.Vascular[,non.op.vars,with = F])[sapply(dt.Vascular[,non.op.vars,with = F], is.integer)]
vars_int <- vars_int[!(vars_int %in% vars_date)]
vars_doub <- names(dt.Vascular[,non.op.vars,with = F])[sapply(dt.Vascular[,non.op.vars,with = F], is.double)]


dt.baseline <- Reduce(function(x, y){merge(x, y, by = "patient_id", all = T)}, 
                      list(dcast(melt(dt.Vascular, id.vars = 'patient_id',measure.vars = vars_char)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Vascular, id.vars = 'patient_id',measure.vars = vars_date)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable, drop = c(T,F)),
                           dcast(melt(dt.Vascular, id.vars = 'patient_id',measure.vars = vars_int)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable),
                           dcast(melt(dt.Vascular, id.vars = 'patient_id',measure.vars = vars_doub)[!is.na(value), .SD[1], by = .(patient_id, variable)], patient_id ~ variable)))

dt.Vascular <- dt.baseline[aggregate_operations, on = "patient_id"]

#####
dt.update <-rbindlist(list(dt.Abdo,
                           dt.Cardiac,
                           dt.Obstetrics,
                           dt.Orthopaedic,
                           dt.Thoracic,
                           dt.Vascular), fill=TRUE) 
rm(dt.Abdo,
   dt.Cardiac,
   dt.Obstetrics,
   dt.Orthopaedic,
   dt.Thoracic,
   dt.Vascular)

dt.update[,dateofbirth := (data.table::as.IDate(paste0(dob,'-15')))]
dt.update[dereg_date != "",gp.end := data.table::as.IDate(paste0(dereg_date,'-15'))]
dt.update[, imd := as.numeric(imd)]
dt.update[, imd5 := cut(imd, breaks = seq(-1,33000,33000/5),  include.lowest = T, ordered_result = F)]


data.table::setkey(dt.update,patient_id)


dt.update[, keep := F]
dt.update[!is.na(dereg_date), dereg_date := data.table::as.IDate(paste0(dereg_date,"-30"), format = "%Y-%m-%d")]
for (x in paste0(procedures,"_date_admitted")) { dt.update[is.na(dereg_date) | x <= dereg_date, keep := T] }
dt.update <- dt.update[keep == T,]

rm(aggregate_operations)

all.days <- names(dt.update)[sapply(dt.update, function(x) inherits(x, 'Date') )]
dt.update[,(all.days) := lapply(.SD, as.IDate),.SDcols = all.days]

data.table::setkey(dt.update,patient_id)

arrow::write_feather(dt.update, sink = here::here("output","update_Sep22.feather"))