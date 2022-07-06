library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
#load(file = here::here("output","cohort_long.RData"))

dt.tv[,`:=`(post.disch.wk1 = discharge.date  ,
            post.disch.wk2 = discharge.date + 7 ,
            post.disch.wk3 = discharge.date + 14 ,
            post.disch.wk4 = discharge.date + 21 ,
            post.disch.wk5 = discharge.date + 28,
            covid.event = event == 1,
            censor.event = event!=1)]

dt.tv.splits <- data.table::melt(dt.tv[is.finite(discharge.date),tail(.SD,1),
                                       .SDcols = c(paste0("post.disch.wk",1:5)),
                                       keyby = c('patient_id')], 
                                 id.vars = c('patient_id'),
                                 value.name = "post.disch.wk",
                                 variable.name = "week.post.disch" )[
                                   order(patient_id, week.post.disch),][,
                                                                        `:=` (tstart = post.disch.wk,
                                                                              week.post.disch = as.numeric(gsub(".*?([0-9]+).*",
                                                                                                                '\\1',
                                                                                                                week.post.disch)))][
                                                                                                                  ,.(patient_id, week.post.disch,tstart)]
dt.tv[,week.post.disch := NA]
dt.tv.splits <- unique(data.table::rbindlist(list(dt.tv.splits,dt.tv[,.(patient_id, week.post.disch, tstart)])))
dt.tv.splits[,post.disch.wk := tstart]
data.table::setkey(dt.tv.splits,patient_id, tstart)
dt.tv[,week.post.disch := NULL]
dt.tv.splits <- dt.tv[dt.tv.splits,,roll = Inf, on = .(patient_id,tstart)]

dates.expand.start.align_(dt = 'dt.tv.splits',
                          start.DTTM = 'tstart',
                          end.DTTM = 'tstop',
                          ID = 'patient_id',
                          merged.DTTM = 'post.disch.wk')

locf.roll_(dt = 'dt.tv.splits',
           ID = 'patient_id',
           start.DTTM = 'tstart',
           group = 'c("patient_id","end.fu")',
           var.cols = 'c("week.post.disch")')

dt.tv.splits[tstart >= admit.date & (!is.finite(discharge.date) | tstop <= discharge.date), week.post.disch := 0 ]

dt.tv.splits[,week.post.disch := as.factor(week.post.disch)]
dt.tv.splits[, los.end := min(los.end, na.rm = T), keyby = .(patient_id, end.fu)]

dt.tv.splits[, `:=`(start = tstart - los.end,
                    end = tstop - los.end)]
dt.tv.splits <- dt.tv.splits[is.finite(los.end) & start>=0 & end <=90]

#Reset outcomes to current end times
dt.tv.splits[,event :=0]
dt.tv.splits[COVIDpositivedate == tstop, event := 1]
dt.tv.splits[emergency_readmitdate  == tstop & event != 1 , event := 2]
dt.tv.splits[date_death_ons == tstop & event != 1, event := 3]

dt.tv.splits[,event.VTE :=0]
dt.tv.splits[post.VTE.date == tstop, event.VTE := 1]
dt.tv.splits[emergency_readmitdate  == tstop & event.VTE != 1, event.VTE := 2]
dt.tv.splits[date_death_ons == tstop & event.VTE != 1, event.VTE := 3]


data.table::setkey(dt.tv.splits, patient_id, end.fu, start)
arrow::write_feather(dt.tv.splits, sink = here::here("output","cohort_postdisch_week_splits.feather"))
#save(dt.tv.splits, file = here::here("output","cohort_postdisch_week_splits.RData"))

