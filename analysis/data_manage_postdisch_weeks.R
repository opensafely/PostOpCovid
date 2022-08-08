library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
#load(file = here::here("output","cohort_long.RData"))

dt.tv[,`:=`(post.op.wk1 = study.start  ,
            post.op.wk2 = study.start + 7 ,
            post.op.wk3 = study.start + 14 ,
            post.op.wk4 = study.start + 21 ,
            post.op.wk5 = study.start + 28,
            covid.event = event == 1,
            censor.event = event!=1)]

# Identify unique date times for week boundaries per patient
dt.tv.splits <- data.table::melt(dt.tv[is.finite(end.fu),tail(.SD,1),
                                       .SDcols = c(paste0("post.op.wk",1:5)),
                                       keyby = c('patient_id','end.fu')], 
                                 id.vars = c('patient_id','end.fu'),
                                 value.name = "post.op.wk",
                                 variable.name = "week.post.op" )[
                                   order(patient_id,end.fu, week.post.op),][,
                                                                        `:=` (tstart = post.op.wk,
                                                                              week.post.op = as.numeric(gsub(".*?([0-9]+).*",
                                                                                                                '\\1',
                                                                                                             week.post.op)))][
                                                                                                                  ,.(patient_id, end.fu, week.post.op,tstart)]
dt.tv[,week.post.op := NA]
dt.tv.splits <- unique(data.table::rbindlist(list(dt.tv.splits,dt.tv[,.(patient_id, end.fu, week.post.op, tstart)])))
dt.tv.splits[,post.op.wk := tstart]
data.table::setkey(dt.tv.splits,patient_id, end.fu, tstart)
dt.tv[,week.post.op := NULL]
dt.tv.splits <- dt.tv[dt.tv.splits,,roll = Inf, on = .(patient_id, end.fu,tstart)]
rm(dt.tv)
gc()
dates.expand.start.align_(dt = 'dt.tv.splits',
                          start.DTTM = 'tstart',
                          end.DTTM = 'tstop',
                          ID = 'patient_id',
                          merged.DTTM = 'post.op.wk')

dates.expand.end.align_(dt = 'dt.tv.splits',
                          start.DTTM = 'tstart',
                          end.DTTM = 'tstop',
                          ID = 'patient_id',
                          merged.DTTM = 'post.op.wk')

locf.roll_(dt = 'dt.tv.splits',
           ID = 'patient_id',
           start.DTTM = 'tstart',
           group = 'c("patient_id","end.fu")',
           var.cols = 'c("week.post.op")')

dt.tv.splits[is.na(week.post.op) & tstart >= admit.date & (!is.finite(discharge.date) | tstop <= discharge.date), week.post.op := 0 ]


dt.tv.splits[,week.post.op := as.factor(week.post.op)]
dt.tv.splits[, los.end := min(los.end, na.rm = T), keyby = .(patient_id, end.fu)]

###################### Post discharge
dt.tv.splits[,`:=`(post.disch.wk1 = discharge.date  ,
            post.disch.wk2 = discharge.date + 7 ,
            post.disch.wk3 = discharge.date + 14 ,
            post.disch.wk4 = discharge.date + 21 ,
            post.disch.wk5 = discharge.date + 28,
            covid.event = event == 1,
            censor.event = event!=1)]

# Identify unique date times for week boundaries per patient
dt.tv.splits2 <- data.table::melt(dt.tv.splits[is.finite(discharge.date),tail(.SD,1),
                                       .SDcols = c(paste0("post.disch.wk",1:5)),
                                       keyby = c('patient_id', 'end.fu')], 
                                 id.vars = c('patient_id', 'end.fu'),
                                 value.name = "post.disch.wk",
                                 variable.name = "week.post.disch" )[
                                   order(patient_id, end.fu, week.post.disch),][,
                                                                        `:=` (tstart = post.disch.wk,
                                                                              week.post.disch = as.numeric(gsub(".*?([0-9]+).*",
                                                                                                                '\\1',
                                                                                                                week.post.disch)))][
                                                                                                                  ,.(patient_id, end.fu, week.post.disch,tstart)]
dt.tv.splits[,week.post.disch := NA]
dt.tv.splits2 <- unique(data.table::rbindlist(list(dt.tv.splits2,dt.tv.splits[,.(patient_id, end.fu, week.post.disch, tstart)])))
dt.tv.splits2[,post.disch.wk := tstart]
data.table::setkey(dt.tv.splits2,patient_id, end.fu, tstart)
dt.tv.splits[,week.post.disch := NULL]
dt.tv.splits <- dt.tv.splits[dt.tv.splits2,,roll = Inf, on = .(patient_id, end.fu,tstart)]
rm(dt.tv.splits2)
dates.expand.start.align_(dt = 'dt.tv.splits',
                          start.DTTM = 'tstart',
                          end.DTTM = 'tstop',
                          ID = 'patient_id',
                          merged.DTTM = 'post.disch.wk')

locf.roll_(dt = 'dt.tv.splits',
           ID = 'patient_id',
           start.DTTM = 'tstart',
           group = 'c("patient_id","end.fu")',
           var.cols = 'c("week.post.disch","postcovid")')

dt.tv.splits[tstart >= admit.date & (!is.finite(discharge.date) | tstop <= discharge.date), week.post.disch := 0 ]

dt.tv.splits[,week.post.disch := as.factor(week.post.disch)]
dt.tv.splits[, los.end := min(los.end, na.rm = T), keyby = .(patient_id, end.fu)]

#Reset outcomes and timings to expanded date times
dt.tv.splits[, `:=`(start = tstart - study.start,
             end = tstop - study.start)]


dt.tv.splits[,discharge.date.locf:= discharge.date]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv.splits',
   min.var.name = 'readmit.end',
   aggregate.cols = 'emergency_readmitdate',
   id.vars = c("patient_id","end.fu"))

dt.tv.splits[, `:=`(start.readmit = tstart - discharge.date.locf,
             end.readmit = tstop - discharge.date.locf)]


## post op COVID cohort----
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[,final.date := covid.end]
dt.tv.splits[is.finite(readmit.end) & readmit.end < final.date & readmit.end > study.start & COVIDreadmission == F, final.date := readmit.end]
dt.tv.splits[is.finite(end.fu) & (end.fu < final.date | !is.finite(final.date)), final.date := end.fu]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv.splits',min.var.name = 'final.date',aggregate.cols = 'final.date',id.vars = c("patient_id","end.fu"))

dt.tv.splits[, postop.covid.cohort := start>=0 & tstop <= final.date & end <= 90]

dt.tv.splits[(postop.covid.cohort) & start ==0  & is.finite(admit.date),any.op.COVID := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv.splits[is.na(any.op.COVID), any.op.COVID := F]
dt.tv.splits[, any.op.COVID := any.op.COVID > 0]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[, any.op.COVID := cummax(any.op.COVID), keyby = .(patient_id, end.fu)]

dt.tv.splits[, postop.covid.cohort := postop.covid.cohort & any.op.COVID == T]

dt.tv.splits[,event :=0]
dt.tv.splits[COVIDpositivedate == tstop & (postop.covid.cohort), event := 1] 
dt.tv.splits[emergency_readmitdate  == tstop & 
                     event != 1 & 
                     readmit.end > study.start & 
                     COVIDreadmission == F &
                     (postop.covid.cohort), event := 2] # Non covid emergency readmission needs to be after admission date (as cannot tell same day admission sequence)
dt.tv.splits[date_death_ons == tstop & event != 1 & (postop.covid.cohort), event := 3]


#### Post discharge events

# VTE events
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[, VTE.end := post.VTE.date]
dt.tv.splits[!is.finite(VTE.end) & (end.fu < post.VTE.date | !is.finite(VTE.end)) , VTE.end := end.fu]

dt.tv.splits[,final.date.VTE := VTE.end]
dt.tv.splits[is.finite(readmit.end) & readmit.end < final.date.VTE & COVIDreadmission == F & readmit.end > study.start, final.date.VTE := readmit.end]
dt.tv.splits[is.finite(end.fu) & (end.fu < final.date.VTE | !is.finite(final.date.VTE)) , final.date.VTE := end.fu]

data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv.splits',min.var.name = 'final.date.VTE',aggregate.cols = 'final.date.VTE',id.vars = c("patient_id","end.fu"))

data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[, postcovid.VTE.cohort := start>=0 & tstop <= final.date.VTE]
dt.tv.splits[(postcovid.VTE.cohort) & start ==0  & is.finite(admit.date),any.op.VTE := rowSums(.SD,na.rm =T)  > 0, .SDcols = c(procedures)]
dt.tv.splits[is.na(any.op.VTE), any.op.VTE := F]
dt.tv.splits[, any.op.VTE := any.op.VTE > 0]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[, any.op.VTE := cummax(any.op.VTE), keyby = .(patient_id, end.fu)]

dt.tv.splits[, postcovid.VTE.cohort := start.readmit > 0 & tstop <= final.date.VTE & any.op.VTE == T] #Date must be after operation date, likely to mean after discharge date as operation date is admit date

dt.tv.splits[,event.VTE :=0]
dt.tv.splits[post.VTE.date == tstop & (postcovid.VTE.cohort), event.VTE := 1]
dt.tv.splits[emergency_readmitdate  == tstop & event.VTE != 1 & COVIDreadmission == F & is.finite(readmit.end) & readmit.end > study.start & (postcovid.VTE.cohort), event.VTE := 2]
dt.tv.splits[date_death_ons == tstop & event.VTE != 1 & (postcovid.VTE.cohort), event.VTE := 3]


## Readmission cohort ----

data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[COVIDreadmission == F & start.readmit > 0 & is.finite(readmit.end) & readmit.end > study.start,final.date.readmit := readmit.end] # Readmission defined as after discharge date in study definition, but if readmitted same day as operation we do not have times to ensure sequence is correct
dt.tv.splits[is.finite(end.fu) & (end.fu < final.date.readmit | !is.finite(final.date.readmit) ) , final.date.readmit := end.fu]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)

min.grp.col_(dt = 'dt.tv.splits',min.var.name = 'final.date.readmit',aggregate.cols = 'final.date.readmit',id.vars = c("patient_id","end.fu"))

dt.tv.splits[, postop.readmit.cohort := start.readmit >= 0 & tstop <= final.date.readmit & end.readmit <= 90]

dt.tv.splits[(postop.readmit.cohort) & start.readmit ==0 ,any.op.readmit := rowSums(.SD,na.rm =T) > 0, .SDcols = c(procedures)]
dt.tv.splits[is.na(any.op.readmit), any.op.readmit := F]
dt.tv.splits[, any.op.readmit := any.op.readmit > 0]
data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)
dt.tv.splits[, any.op.readmit := cummax(any.op.readmit), keyby = .(patient_id, end.fu)]

dt.tv.splits[, postop.readmit.cohort :=  postop.readmit.cohort & any.op.readmit == T]


dt.tv.splits[,event.readmit :=0]
dt.tv.splits[(postop.readmit.cohort) & 
      emergency_readmitdate  == tstop & 
      COVIDpositivedate != tstop &  
      COVIDreadmission == F & 
      readmit.end > study.start & is.finite(readmit.end), event.readmit := 1] # COVID isn't a competing risk its an exposure for non covid readmission
dt.tv.splits[(postop.readmit.cohort) & date_death_ons == tstop & event.readmit != 1, event.readmit := 2]

data.table::setkey(dt.tv.splits,patient_id,tstart,tstop)


data.table::setkey(dt.tv.splits, patient_id, end.fu, start)
arrow::write_feather(dt.tv.splits, sink = here::here("output","cohort_postdisch_week_splits.feather"))
#save(dt.tv.splits, file = here::here("output","cohort_postdisch_week_splits.RData"))

