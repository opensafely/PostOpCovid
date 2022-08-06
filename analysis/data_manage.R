# detach("package:here", unload = TRUE)
# setwd("C:\\Users\\mczcjc\\Documents\\GitHub\\PostOpCovid")
# library(here)
# detach("package:here", unload = TRUE)
# setwd("P:\\GitHub\\PostOpCovid")
# library(here)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)
source(here::here("analysis","Utils.R"))

index_date <- data.table::as.IDate("2020-02-01")

dt <- data.table::fread(here::here("output", "input.csv"))
dt.COD <- data.table::fread(here::here("output", "input_COD.csv"))

data.table::setkey(dt,patient_id)
data.table::setkey(dt.COD,patient_id)

dt <- dt.COD[,.(patient_id, death_underlying_cause_ons)][dt,]

####
# Basic counts and descriptions----
####
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

dt[,dateofbirth := (data.table::as.IDate(paste0(dob,'-15')))]
dt[dereg_date != "",gp.end := data.table::as.IDate(paste0(dereg_date,'-15'))]
dt[, imd := as.numeric(imd)]
dt[, imd5 := cut(imd, breaks = seq(-1,33000,33000/5),  include.lowest = T, ordered_result = F)]

####
# Multiple operations per row. Reshape to long format ----
####

## Clean invalid admission sequences ----
op.admit.vars <- sort(paste0(outer(procedures,paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures,paste0("_",1:5),paste0),'_date_discharged'))
for (i in 1:length(op.admit.vars)) {
  dt[get(op.admit.vars[i])>get(op.discharge.vars[i]), 
     (names(dt)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                            pattern = '_date_admitted',
                                            replacement = ""),'*'),
                      x = names(dt))]) := NA]
}
## Reshape repeated procedures within a specialty ----
repeated.vars <- names(dt)[grepl(pattern = paste(procedures,collapse = '|'),x = names(dt))]
non.op.vars <- names(dt)[!grepl(pattern = paste(procedures,collapse = '|'),x = names(dt))]

aggregate_operations <- data.table::melt(dt, id.var = 'patient_id',
                                         measure = patterns(as.vector(outer(paste0(procedures,'_[0-9]'),
                                                                            unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = paste0(repeated.vars,'$'))),
                                                                            paste0))),
                                         variable.name = 'op.number',
                                         value.name = as.vector(outer(procedures,unique(gsub(pattern = paste0(as.vector(outer(procedures,paste0("_",1:5),paste0)),collapse = "|"),replacement = "",x = repeated.vars)),paste0)))[order(patient_id,op.number)]
dt <- dt[,non.op.vars, with = F][aggregate_operations, on = "patient_id"]
rm(aggregate_operations)

dt[, keep := F]
dt[!is.na(dereg_date), dereg_date := data.table::as.IDate(paste0(dereg_date,"-30"), format = "%Y-%m-%d")]
for (x in paste0(procedures,"_date_admitted")) { dt[is.na(dereg_date) | x <= dereg_date, keep := T] }
dt <- dt[keep == T,]
lapply(paste0(procedures,"_date_admitted"), function(x) dt[is.finite(get(x)),.N]) # Check numbers in log

## Waves defined from ONS reports: ----
#https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveycharacteristicsofpeopletestingpositiveforcovid19uk/latest
#https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveyantibodyandvaccinationdatafortheuk/1june2022
dt[,(paste("admit.wave.",procedures, sep ="")) := lapply(.SD, function(x) cut(as.numeric(x), breaks = c(as.numeric(data.table::as.IDate("2020-01-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2020-12-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-05-17", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-12-19", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate(Sys.Date(), format = "%Y-%m-%d"))),
                                                                              labels = c("Wave_1","Wave_2","Wave_3","Wave_4"),
                                                                              ordered = T)), 
                                                         .SDcols = c(paste0(procedures,"_date_admitted"))]

## Table to check numbers before reshaping for debugging only- not run now
#dt[, (paste0(x,"post.VTE")) := ((!is.na(.SD[,3]) &  
#                        .SD[,3] <= .SD[,1] + 90 & .SD[,3] >= .SD[,1]) | 
#                       ((!is.na(.SD[,4]) &  
#                           .SD[,4] <= .SD[,1] + 90 & .SD[,4] >= .SD[,1]))) & 
#        (!is.na(.SD[,5]) &  
#           .SD[,5] <= .SD[,1] + 90 & .SD[,5] >= .SD[,1]), 
#   .SDcols = paste0(x,c("_date_admitted","_date_discharged","_VTE_GP_date","_VTE_HES_date_admitted","_anticoagulation_prescriptions_date"))]  # events flagged at end of episode
#}

# demo.waves.tab <- lapply(procedures, function(proc) { 
#   t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  {
#     cbind(dt[is.finite(get(paste0(proc,'_date_admitted'))) & get(paste0('admit.wave.',proc)) == x,.("Procedures" = .N,
#                                  "Patients" = length(unique(patient_id)),
#                                  "Male" = round(mean(sex=='M'),digits = 2),
#                                  "Age (IQR)" = paste(round(quantile(as.numeric(get(paste0(proc,'_date_admitted')) - as.numeric(dateofbirth))/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
#                                  "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
#                                  "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
#                                  "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_1 >= 14),digits = 2),
#                                  "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_2 >= 14),digits = 2),
#                                  "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_3 >= 14),digits = 2),
#                                  "Current Cancer"  = round(mean(substr(get(paste0(proc,'_primary_diagnosis')),1,1) =='C'), digits = 2),
#                                  "Emergency" = round(mean(substr(get(paste0(proc,'_admission_method')),1,1) == "2"),digits = 2),
#                                  "Length of stay (IQR)" =  paste(round(quantile((get(paste0(proc,'_date_discharged')) - get(paste0(proc,'_date_admitted'))),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
#                                  "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - get(paste0(proc,'_date_admitted')) <= 90),digits = 2),
#                                  "90 day COVID-19" = round(mean(is.finite(get(paste0(proc,'_date'))) & get(paste0(proc,'_date')) - get(paste0(proc,'_date_admitted')) <= 90  & get(paste0(proc,'_date')) - get(paste0(proc,'_date_admitted')) >=0),digits = 2),
#                                  "90 day VTE" = round(mean(get(paste0(proc,'post.VTE')), na.rm = T),digits = 2))],
#                                t(dt[is.finite(get(paste0(proc,'_date_admitted'))) & get(paste0('admit.wave.',proc)) == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))]),
#                                t(dt[get(paste0('admit.wave.',proc)) == x,.N, by = .(.grp = eval(parse(text = paste0(proc,'_primary_diagnosis'))))][order(-N), do.call(paste,c(.SD, sep = ": "))][1:5])                              
#                                )})), fill = T))})
# 
# demo.waves.tab
# 
# for(i in 1:length(procedures)) print(xtable::xtable(demo.waves.tab[[i]]), type = 'html', here::here("output",paste0("table1",procedures[i],".html")))
# rm(demo.waves.tab)

####
# Reshape data to long time varying cohort per procedure ----
####

## Baseline fixed factors ----
fixed <- c('patient_id','dob','sex','bmi' ,'region', 'imd','date_death_ons')

## Exposures unrelated to operation times ----
time.cols <- c(paste0("covid_vaccine_dates_",1:3),c(names(dt)[grep("^pre",names(dt))])) 

## Time varying exposures associated with each procedure date:  need to start from beginning of row time period ----

proc.tval.stubs <- c('_admission_method',
                     '_primary_diagnosis',
                     '_days_in_critical_care',
                     '_case_category',
                     '_recent_case_category',
                     '_previous_case_category',
                     '_HipReplacement_HES_binary_flag',
                     '_KneeReplacement_HES_binary_flag',
                     '_Cholecystectomy_HES_binary_flag',
                     '_Colectomy_HES_binary_flag')
proc.tval.cols <- paste(rep(procedures,each = length(proc.tval.stubs)),proc.tval.stubs, sep ="")

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
dt[,(c(proc.time.cols.start,proc.time.cols.end)) := lapply(.SD,as.numeric), .SDcols = c(proc.time.cols.start,proc.time.cols.end)]
dt[,(time.cols) := lapply(.SD,as.numeric), .SDcols = time.cols]
dt[,date_death_ons := as.numeric(date_death_ons)]
dt[,gp.end := as.numeric(gp.end)]

# Variable for earliest of 90 days post discharge or death for maximum end of follow up time
## 90 days Post discharge ----
dt[,(paste0(procedures,"_end_fu")) := lapply(.SD, 
                                             function(x) data.table::fifelse(is.finite(date_death_ons) & x+90 > date_death_ons, 
                                                                             date_death_ons,
                                                                             x+90)),
   .SDcols = paste0(procedures, '_date_discharged')]

dt[,(paste0(procedures,"_end_fu")) := lapply(.SD,
                                             function(x) data.table::fifelse(is.finite(gp.end) & x > gp.end,
                                                                             gp.end,
                                                                             x)),
   .SDcols = paste0(procedures, '_end_fu')]

## 90 days Post operative to ensure time break to allow censoring here  ----
dt[,(paste0(procedures,"_end_fu90")) := lapply(.SD, 
                                               function(x) data.table::fifelse(is.finite(date_death_ons) & x+90 > date_death_ons, 
                                                                               date_death_ons,
                                                                               x+90)),
   .SDcols = paste0(procedures, '_date_admitted')]

dt[,(paste0(procedures,"_end_fu90")) := lapply(.SD, 
                                               function(x) data.table::fifelse(is.finite(gp.end) & x > gp.end, 
                                                                               gp.end,
                                                                               x)),
   .SDcols = paste0(procedures, '_end_fu90')]

## 30 days Post operative to ensure time break to allow censoring here ----
dt[,(paste0(procedures,"_end_fu30")) := lapply(.SD, 
                                               function(x) data.table::fifelse(is.finite(date_death_ons) & x+30 > date_death_ons,
                                                                               date_death_ons,
                                                                               x+30)),
   .SDcols = paste0(procedures, '_date_admitted')]

dt[,(paste0(procedures,"_end_fu30")) := lapply(.SD,
                                               function(x) data.table::fifelse(is.finite(gp.end) & x > gp.end,
                                                                               gp.end,
                                                                               x)),
   .SDcols = paste0(procedures, '_end_fu30')]


## Final date for per patient taking into account final patient in GP database and currrent study end date ----

max.grp.col_(dt = 'dt', max.var.name = 'max.date', aggregate.cols = paste0(procedures,"_end_fu"), id.vars = 'patient_id') 
min.grp.col_(dt = 'dt', min.var.name = 'min.date', aggregate.cols = c(time.cols,proc.time.cols.start,proc.time.cols.end), id.vars = 'patient_id') 


dt[is.finite(gp.end) & max.date > gp.end, max.date := gp.end]
dt[!is.finite(max.date), max.date := as.numeric(data.table::as.IDate('2022-02-01'))]
dt[,tstart := do.call(pmin, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_admitted")]
dt[op.number == 1 ,tstart:= min.date]


## All times for long cohort table ----
dt.dates.wide <- dt[,.SD, .SDcols = c('patient_id',
                                 time.cols,
                                 'gp.end',
                                 proc.time.cols.start, 
                                 proc.time.cols.end,
                                 paste(procedures,"_end_fu",sep = ""),
                                 paste(procedures,"_end_fu30",sep = ""),
                                 paste(procedures,"_end_fu90",sep = ""))]
dt.dates.wide[, gp.end := as.numeric(gp.end)]


# Roll dates into main data ----

# melt dates to long
dt.dates.long <- unique(data.table::melt(dt.dates.wide, id.vars = 'patient_id',na.rm = T, value.name = 'tstart')[is.finite(tstart),c(1,3)])
data.table::setkey(data.table::setDT(dt),patient_id, tstart)
data.table::setkey(dt.dates.long, patient_id, tstart)

# Roll dates into data with locf for all variables, (and nocb for any earlier times)
dt.tv <- dt[dt.dates.long,,rollends = c(T,T), roll = Inf, on = c('patient_id','tstart'), mult = 'all']
rm(dt)
data.table::setkey(dt.tv,patient_id, tstart)

# Set date end for each period
dt.tv[,tstop := c(tstart[-1],NA)] 
dt.tv[c(F,patient_id[c(-1,-.N)] != patient_id[c(-1,-2)],T), tstop:=max.date]

####
# Align time index across all records within patient and procedure----
####

## Drop dates from copied rows 

## tstart aligned dates ----
data.table::setkey(dt.tv,patient_id,tstart,tstop)

for (v in c(proc.time.cols.start, time.cols))  dt.tv[get(v) != tstart, (v) := NA]

for (i in 1:length(proc.tval.cols))  dt.tv[get(paste0(strsplit(proc.tval.cols[i], "[_]")[[1]][1],"_date_admitted")) != tstart, (proc.tval.cols[i]) := NA]

## tstop aligned dates ----
for (i in 1:length(procedures))  dt.tv[get(paste0(procedures[i],
                                                  "_emergency_readmit_date_admitted")) != tstop , 
                                       (c(paste0(procedures[i],"_emergency_readmit_date_admitted"),
                                          paste0(procedures[i],"_emergency_readmit_primary_diagnosis"))) := NA]


for (i in c(proc.time.cols.end,
            paste0(procedures,"_end_fu"),
            paste0(procedures,"_end_fu30"),
            paste0(procedures,"_end_fu90"),
            "gp.end"))  dt.tv[get(v) != tstop, (v) := NA]


# Copy gp end of follow up across all patient time
data.table::setorder(dt.tv,patient_id, tstart, tstop)
max.grp.col_(dt = 'dt.tv', max.var.name = 'gp.end', aggregate.cols = 'gp.end', id.vars = 'patient_id') 

dt.tv[gp.end == 0, gp.end := Inf]

####
# Coalesce shared dates across operation types ----
####

data.table::setkey(dt.tv,patient_id, tstart, tstop)
admission.dates <- c('admit.date','discharge.date','end.fu') # key dates to define

## Admission dates from different procedures to create continuous record, only present when valid ----
dt.tv[,admit.date := do.call(pmin, c(.SD, na.rm = T)),
      .SDcols = paste0(procedures,"_date_admitted")]
dt.tv[!is.finite(admit.date) | admit.date != tstart, admit.date := NA]

## End of follow up from different procedures to create continuous record, only present when valid ----
dt.tv[,end.fu := do.call(pmin, c(.SD, na.rm = T)), 
      .SDcols = c(paste0(procedures,"_end_fu"))] ## gp.end in end_fu already
dt.tv[!is.finite(end.fu) | end.fu!= tstop, end.fu := NA]

dt.tv[,end.fu30 := do.call(pmin, c(.SD, na.rm = T)), 
      .SDcols = c(paste0(procedures,"_end_fu30"))] ## gp.end in end_fu already
dt.tv[!is.finite(end.fu30) | end.fu30!= tstop, end.fu30 := NA]

dt.tv[,end.fu90 := do.call(pmin, c(.SD, na.rm = T)), 
      .SDcols = c(paste0(procedures,"_end_fu90"))] ## gp.end in end_fu already
dt.tv[!is.finite(end.fu90) | end.fu90!= tstop, end.fu90 := NA]

dt.tv[,discharge.date := do.call(pmax, c(.SD, na.rm = T)), 
      .SDcols = paste0(procedures,"_date_discharged")]
dt.tv[!is.finite(discharge.date) | discharge.date!=tstop, discharge.date := NA]
data.table::setkey(dt.tv,patient_id,tstart,tstop)

####
## Operation spell periods ----
####


## Next admission period / procedure censor end of each spell
data.table::setkey(dt.tv,patient_id, tstart, tstop)
dt.tv[,next.admit.date := data.table::shift(admit.date, n = 1L, type =  'lead'), by = .(patient_id)]
nocb.roll_(dt = 'dt.tv', ID = 'patient_id', 
           start.DTTM = 'tstart', group = 'c("patient_id")',
           var.cols = 'c("next.admit.date")')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[next.admit.date < end.fu, end.fu := next.admit.date]
dt.tv[next.admit.date < end.fu30, end.fu30 := next.admit.date]
dt.tv[next.admit.date < end.fu90, end.fu90 := next.admit.date]

## Roll end of period dates backwards to end of previous episodes to define each post procedure period ----
data.table::setkey(dt.tv,patient_id, tstart, tstop)
nocb.roll_(dt = 'dt.tv', ID = 'patient_id', 
           start.DTTM = 'tstart', group = 'c("patient_id")',
           var.cols = 'c("discharge.date","end.fu")')

dt.tv[discharge.date > end.fu, discharge.date := NA]
dt.tv[tstart >= discharge.date, discharge.date := NA]

## Start of follow up (each enter study) per patient and procedure
data.table::setkey(dt.tv,patient_id, tstart, tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'study.start',
             aggregate.cols = 'admit.date',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(study.start), study.start := NA]

## Roll start of procedure periods forward to define beginning of each post procedure period. Include discharge dates and end of fu dates to work out data to drop later ----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv', ID = 'patient_id', 
           start.DTTM = 'tstart', group = 'c("patient_id")',
           var.cols = 'c(admission.dates)')


## Remove admit & discharge dates from outside admission period. Study start and endfu define total exposure for each observation period.
dt.tv[admit.date > discharge.date | is.na(admit.date) | admit.date > tstart | discharge.date < tstop, c('admit.date','discharge.date') := NA]
dt.tv[is.na(admit.date) , (as.vector(outer(procedures,c('_admission_method','_primary_diagnosis',
                                                        '_days_in_critical_care'), paste0)))  := NA]
dt.tv <- dt.tv[tstop <= end.fu,]

# Pre study start not dropped yet as contains pre operative exposure information

# Set infinite to missing to avoid invalid comparisons - not a problem with rolls /sorting min & max
dt.tv[!is.finite(admit.date), admit.date := NA]
dt.tv[!is.finite(end.fu), end.fu := NA]
dt.tv[!is.finite(end.fu30), end.fu30 := NA]
dt.tv[!is.finite(end.fu90), end.fu90 := NA]
dt.tv[!is.finite(discharge.date), discharge.date := NA]
dt.tv[!is.finite(study.start), study.start := NA]

####
# Defining operation for each post op period--------------
####

for (i in 1:length(procedures)) dt.tv[, (procedures[i]) := is.finite(get(paste0(procedures[i],"_date_admitted"))) &
                                        is.finite(admit.date) &
                                        admit.date <= get(paste0(procedures[i],"_date_admitted")) & 
                                         (is.finite(end.fu) | 
                                           end.fu >= get(paste0(procedures[i],"_date_admitted")))]

dt.tv[,(procedures) := lapply(.SD, function(x) data.table::fifelse(x==0,NA,x)), .SDcols = c(procedures)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv', ID = 'patient_id', start.DTTM = 'tstart', group = 'c("patient_id","end.fu")', var.cols = 'c(procedures)')
for (i in 1:length(procedures)) dt.tv[!is.finite(get(procedures[i])), (procedures[i]) := F]

## Other time varying information coalesced for each procedure spell----

data.table::setkey(dt.tv,patient_id,tstart,tstop)

for(v in c(proc.tval.stubs,'_emergency_readmit_primary_diagnosis')) { 
  data.table::setkey(dt.tv,patient_id,tstart,tstop)
  max.grp.col_(dt = 'dt.tv',
               max.var.name = gsub("*_HES_binary_flag$","",gsub("^_*","",v)),
               aggregate.cols = paste0(procedures,v),
               id.vars = c("patient_id","end.fu"))
  dt.tv[,(paste0(procedures,v)) := NULL]
}


for(v in c(proc.time.stubs.start[!(proc.time.stubs.start == '_date_admitted')],
           proc.time.stubs.end[!(proc.time.stubs.end %in% c('_date_discharged','_end_fu'))])) { 
  data.table::setkey(dt.tv,patient_id,tstart,tstop)
  min.grp.col_(dt = 'dt.tv',
               min.var.name = gsub("^_*","",v),
               aggregate.cols = paste0(procedures,v),
               id.vars = c("patient_id","end.fu"))
}


dt.tv[,Current.Cancer := substr(primary_diagnosis,1,1) =='C']
dt.tv[is.na(Current.Cancer), Current.Cancer := F]

dt.tv[,(c(proc.time.cols.start,
          proc.time.cols.end, 
          proc.tval.cols,
          paste0(procedures,'_emergency_readmit_primary_diagnosis'),
          paste0(procedures,'_end_fu'),
          paste0(procedures,'post_VTE'))) := NULL]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv', ID = 'patient_id', start.DTTM = 'tstart', group = 'c("patient_id","end.fu")', var.cols = 'c(time.cols)')

####
# Pre operative risk factors----
####
dt.tv[,age := floor((tstart - as.numeric(as.Date(paste0(dob,'-15'))))/365.25)]
max.age <- max(dt.tv$age,na.rm = T)
dt.tv[,age.cat := cut(age, breaks = c(1,50,70,80,max.age),ordered_result = F , right = T, include.lowest = T)]

dt.tv[, imd5 := cut(imd, breaks = seq(0,33000,33000/5),  include.lowest = T, ordered_result = F)]
if(sum(is.na(dt.tv$imd5))!=0) {
  levels(dt.tv$imd5) <- c(levels(dt.tv$imd5),"Missing")
  dt.tv[is.na(imd5) , imd5 := "Missing"]
}

dt.tv[, bmi.cat := cut(bmi, breaks = c(0,18,24,29,100),  include.lowest = T, ordered_result = F)]
if(sum(is.na(dt.tv$bmi.cat))!=0) {
  levels(dt.tv$bmi.cat) <- c(levels(dt.tv$bmi.cat),"Missing")
  dt.tv[is.na(bmi.cat) , bmi.cat := "Missing"]
}

dt.tv[, region := as.factor(region)]
if(sum(is.na(dt.tv$region))!=0) {
  levels(dt.tv$region) <- c(levels(dt.tv$region),"Missing")
  dt.tv[is.na(region) , region := "Missing"]
}

## Charlson index at time of operation----

data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv[, Charlson := is.finite(pre_MI_GP) +
        is.finite(pre_CCF_GP) +
        is.finite(pre_PVD_GP) +
        is.finite(pre_Stroke_GP) +
        is.finite(pre_Dementia_GP) +
        is.finite(pre_Respiratory_GP) +
        is.finite(pre_RA_SLE_Psoriasis_GP) +
        is.finite(pre_Ulcer_or_bleed_GP) +
        is.finite(pre_all_liver_GP) + 
        is.finite(pre_Cirrhosis_GP)*2 + # counted in all_liver_GP too
        is.finite(pre_all_diabetes_GP) +
        is.finite(pre_Diabetic_Complications_GP) + # counted in diabetes too
        is.finite(pre_Other_Neurology_GP)*2 +
        (is.finite(pre_CKD_3_5_GP) | is.finite(pre_Renal_GP))*2 +
        is.finite(pre_Non_Haematology_malignancy_GP)*2 +
        is.finite(pre_Haematology_malignancy_GP)*2 +
        is.finite(pre_Metastases_GP)*6 +
        is.finite(pre_HIV_GP)*6]  

dt.tv[,Charl12 := cut(Charlson, breaks = c(0,1,2,100),  include.lowest = T, labels = c("None","Single","Multiple or Severe"), ordered = F)]
if(sum(is.na(dt.tv$Charl12))!=0) {
  levels(dt.tv$Charl12) <- c(levels(dt.tv$Charl12),"Missing")
  dt.tv[is.na(Charl12) , Charl12 := "Missing"]
}

## Elective or emergency operations----
dt.tv[,Emergency := substr(as.character(admission_method),1,1)=="2"]
dt.tv[is.na(Emergency), Emergency := F]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'Emergency',aggregate.cols = 'Emergency',id.vars = c("patient_id","end.fu"))

## Vaccination status - 14 days post date as effective----
dt.tv[, vaccination.status := is.finite(covid_vaccine_dates_1) + is.finite(covid_vaccine_dates_2) + is.finite(covid_vaccine_dates_3)]
dt.tv[,vaccination.status.factor := factor(vaccination.status,  ordered = F)]
dt.tv[is.na(vaccination.status.factor), vaccination.status.factor := 0]

####
# Post operative outcomes----
####

## Length of stay----
dt.tv[, LOS := discharge.date - admit.date]


data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = "LOS",
             aggregate.cols = 'LOS',
             id.vars = c("patient_id","end.fu"))
dt.tv[,LOS.bin := LOS > 7]

dt.tv[,discharged := is.finite(discharge.date) & discharge.date == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'los.end',aggregate.cols = 'discharge.date',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(los.end), los.end := end.fu]

## VTE----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, anticoagulation_prescriptions_date := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = 'anticoagulation_prescriptions_date']
dt.tv[, anticoagulation_prescriptions_date := lapply(.SD, data.table::nafill, type = "locf"), by = patient_id, .SDcols = 'anticoagulation_prescriptions_date']

dt.tv[is.finite(death_underlying_cause_ons) & death_underlying_cause_ons %in% c('O225','O873','G08','I636','I676','I80','I801','I802','I803','O223','O871','I820',
'I26','I260','I269','O082','O882','I81','I823','I808','I809','I82','I821','I828','I829','I822'), VTE_death := date_death_ons]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',
             min.var.name = "post.VTE.date",
             aggregate.cols = c('VTE_GP_date','VTE_HES_date_admitted','date_death_ons'),
             id.vars = c("patient_id","end.fu"))

dt.tv[is.finite(anticoagulation_prescriptions_date) & 
         is.finite(post.VTE.date) &
        (anticoagulation_prescriptions_date < post.VTE.date - 15 |
        anticoagulation_prescriptions_date > post.VTE.date + 90) , post.VTE.date := NA]

dt.tv[, post.VTE := is.finite(post.VTE.date) & post.VTE.date == tstop]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'post.VTE.date',aggregate.cols = 'post.VTE.date',id.vars = c("patient_id","end.fu"))
dt.tv[,postVTEany := cumsum(post.VTE), by = .(patient_id, end.fu)]

dt.tv[, VTE.end := post.VTE.date]
dt.tv[!is.finite(VTE.end) | post.VTE.date > end.fu, VTE.end := end.fu]

## Covid-19----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
names(dt.tv)[names(dt.tv)=='date'] <- 'COVIDpositivedate'

dt.tv[is.na(COVIDpositivedate) & death_underlying_cause_ons %in% c('U071','U072'), COVIDpositivedate := date_death_ons]

dt.tv[(!is.na(emergency_readmit_primary_diagnosis) & 
emergency_readmit_primary_diagnosis %in% c('U071','U072')) &
 tstop < COVIDpositivedate, COVIDpositivedate := tstop] ## Primary readmission diagnosis was COVID then count at admission time

 data.table::setkey(dt.tv,patient_id,tstart,tstop)
  min.grp.col_(dt = 'dt.tv',
               min.var.name = 'COVIDpositivedate',
               aggregate.cols = 'COVIDpositivedate',
               id.vars = c("patient_id","end.fu"))
dt.tv[,COVIDpositive := is.finite(COVIDpositivedate) & COVIDpositivedate == tstop]
dt.tv[ !is.finite(COVIDpositive), COVIDpositive  := F]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[order(tstart),postcovid := cummax(COVIDpositive), by = .(patient_id, end.fu)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, covid.end := COVIDpositivedate]
dt.tv[!is.finite(covid.end), covid.end := end.fu]

names(dt.tv)[names(dt.tv)=='recent_date'] <- 'recentCOVIDpositivedate'
dt.tv[,recentCOVID := as.numeric(is.finite(recentCOVIDpositivedate)) ]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'recentCOVID',aggregate.cols = 'recentCOVID',id.vars = c("patient_id","end.fu"))
dt.tv[, recentCOVID := recentCOVID==1]
dt.tv[!is.finite(recentCOVID), recentCOVID := F]

names(dt.tv)[names(dt.tv)=='previous_date'] <- 'previousCOVIDpositivedate'
dt.tv[,previousCOVID := as.numeric(is.finite(previousCOVIDpositivedate)) ]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'previousCOVID',aggregate.cols = 'previousCOVID',id.vars = c("patient_id","end.fu"))
dt.tv[, previousCOVID := previousCOVID==1]
dt.tv[!is.finite(previousCOVID), previousCOVID := F]

## Readmissions----
## Ensure readmissions are not mistaken as the initial admission. Discharged date already rolled across end.fu so can do this across all records immediately
dt.tv[emergency_readmit_date_admitted <= discharge.date, `:=`(emergency_readmit_date_admitted = NA,
                                                               emergency_readmit_primary_diagnosis = NA) ]

## Label COVID readmissions if within 2 days of admission in line with Government definition or coded as a primary diagnosis as cannot tell where would occur in readmission
dt.tv[,COVIDreadmission := (!is.na(emergency_readmit_primary_diagnosis) & 
      emergency_readmit_primary_diagnosis %in% c('U071','U072') & emergency_readmit_date_admitted > discharge.date) | 
      (is.finite(COVIDpositivedate) & 
      is.finite(emergency_readmit_date_admitted) & 
      emergency_readmit_date_admitted > discharge.date &
        COVIDpositivedate >= emergency_readmit_date_admitted & 
        COVIDpositivedate <= emergency_readmit_date_admitted + 2)]

dt.tv[!is.finite(COVIDreadmission), COVIDreadmission := F]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
   max.var.name = 'COVIDreadmission',
   aggregate.cols = 'COVIDreadmission',
   id.vars = c("patient_id","end.fu"))

dt.tv[COVIDreadmission == T, `:=`(emergency_readmit_date_admitted = NA,
                                  emergency_readmit_primary_diagnosis = NA)] ## COVID readmission outcomes are in COVID outcomes, keep readmissions as non covid readmissions

dt.tv[,emergency_readmit  := is.finite(emergency_readmit_date_admitted) & 
        COVIDreadmission == F &
        emergency_readmit_date_admitted == tstop]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',
   min.var.name = 'readmit.end',
   aggregate.cols = 'emergency_readmit_date_admitted',
   id.vars = c("patient_id","end.fu"))

names(dt.tv)[names(dt.tv)=='emergency_readmit_date_admitted'] <- 'emergency_readmitdate'

## Mortality----
## date_death_ons part of definition of end_fu so will be end of final row when in follow up period
dt.tv[,died := is.finite(date_death_ons) & tstop == date_death_ons]
dt.tv[is.na(died), died := 0]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv <- dt.tv[!(tstart < study.start | tstop > end.fu) & !is.na(age.cat),]
dt.tv[, year := data.table::year(data.table::as.IDate(admit.date))]

## Waves Redefine in long table from ONS reports:  ----
#https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveycharacteristicsofpeopletestingpositiveforcovid19uk/latest
#https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveyantibodyandvaccinationdatafortheuk/1june2022

dt.tv[,wave := cut(study.start, breaks = c(as.numeric(data.table::as.IDate("2020-01-01", format = "%Y-%m-%d")), ## Wuhan 
                                           as.numeric(data.table::as.IDate("2020-12-01", format = "%Y-%m-%d")), ## alpha
                                           as.numeric(data.table::as.IDate("2021-05-17", format = "%Y-%m-%d")), ## delta
                                           as.numeric(data.table::as.IDate("2021-12-19", format = "%Y-%m-%d")), ## Omicron
                                           as.numeric(data.table::as.IDate(Sys.Date(), format = "%Y-%m-%d"))),
                   labels = c("Wave_1","Wave_2","Wave_3","Wave_4"),
                   include.lowest = T,
                   right = T,
                   ordered = F)]

# Define cohorts ----
## Restart clock with each procedure----
dt.tv[, `:=`(start = tstart - study.start,
             end = tstop - study.start)]
min.grp.col_(dt = 'dt.tv[start >= 0,]',min.var.name = 'discharge.start',aggregate.cols = 'discharge.date',id.vars = c("patient_id","end.fu"))

data.table::setkey(dt.tv,patient_id,tstart,tstop)

# Check all study episodes have a procedure
dt.tv[start ==0  & is.finite(admit.date),any.op := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op), any.op := F]
dt.tv[, any.op := any.op > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op := cummax(any.op), keyby = .(patient_id, end.fu)]

## post op COVID cohort----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,final.date := covid.end]
dt.tv[is.finite(readmit.end) & readmit.end < final.date & readmit.end > study.start & COVIDreadmission == F, final.date := readmit.end]
dt.tv[is.finite(end.fu) & end.fu < final.date, final.date := end.fu]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date',aggregate.cols = 'final.date',id.vars = c("patient_id","end.fu"))

dt.tv[, postop.covid.cohort := start>=0 & tstop <= final.date & end <= 90]

dt.tv[(postop.covid.cohort) & start ==0  & is.finite(admit.date),any.op.COVID := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op.COVID), any.op.COVID := F]
dt.tv[, any.op.COVID := any.op.COVID > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op.COVID := cummax(any.op.COVID), keyby = .(patient_id, end.fu)]

dt.tv[, postop.covid.cohort := start>=0 & tstop <= final.date & end <= 90 & any.op.COVID == T]

dt.tv[,event :=0]
dt.tv[COVIDpositivedate == tstop & (postop.covid.cohort), event := 1] 
dt.tv[emergency_readmitdate  == tstop & 
                     event != 1 & 
                     readmit.end > study.start & 
                     COVIDreadmission == F &
                     (postop.covid.cohort), event := 2] # Non covid emergency readmission needs to be after admission date (as cannot tell same day admission sequence)
dt.tv[date_death_ons == tstop & event != 1 & (postop.covid.cohort), event := 3]


## VTE cohort post Covid----

dt.tv[,discharge.date.locf:= discharge.date]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

min.grp.col_(dt = 'dt.tv',min.var.name = 'discharge.date.locf',aggregate.cols = 'discharge.date.locf',id.vars = c("patient_id","end.fu"))


dt.tv[, `:=`(start.readmit = tstart - discharge.date.locf,
             end.readmit = tstop - discharge.date.locf)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,final.date.VTE := VTE.end]
dt.tv[is.finite(readmit.end) & readmit.end < final.date.VTE & COVIDreadmission == F & readmit.end > study.start, final.date.VTE := readmit.end]
dt.tv[is.finite(end.fu) & end.fu < final.date.VTE, final.date.VTE := end.fu]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date.VTE',aggregate.cols = 'final.date.VTE',id.vars = c("patient_id","end.fu"))

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, postcovid.VTE.cohort := start>=0 & tstop <= final.date.VTE]
dt.tv[(postcovid.VTE.cohort) & start ==0  & is.finite(admit.date),any.op.VTE := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op.VTE), any.op.VTE := F]
dt.tv[, any.op.VTE := any.op.VTE > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op.VTE := cummax(any.op.VTE), keyby = .(patient_id, end.fu)]

dt.tv[, postcovid.VTE.cohort := start.readmit > 0 & tstop <= final.date.VTE & any.op.VTE == T] #Date must be after operation date, likely to mean after discharge date as operation date is admit date

dt.tv[,event.VTE :=0]
dt.tv[post.VTE.date == tstop & (postcovid.VTE.cohort), event.VTE := 1]
dt.tv[emergency_readmitdate  == tstop & event.VTE != 1 & COVIDreadmission == F & readmit.end > study.start & (postcovid.VTE.cohort), event.VTE := 2]
dt.tv[date_death_ons == tstop & event.VTE != 1 & (postcovid.VTE.cohort), event.VTE := 3]



## Readmission cohort ----

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[COVIDreadmission == F & start.readmit > 0 & readmit.end > study.start,final.date.readmit := readmit.end] # Readmission defined as after discharge date in study definition, but if readmitted same day as operation we do not have times to ensure sequence is correct
dt.tv[is.finite(end.fu) & end.fu < final.date.readmit , final.date.readmit := end.fu]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date.readmit',aggregate.cols = 'final.date.readmit',id.vars = c("patient_id","end.fu"))

dt.tv[, postop.readmit.cohort := start.readmit >= 0 & tstop <= final.date.readmit & end.readmit <= 90]

dt.tv[(postop.readmit.cohort) & start.readmit ==0 ,any.op.readmit := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op.readmit), any.op.readmit := F]
dt.tv[, any.op.readmit := any.op.readmit > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op.readmit := cummax(any.op.readmit), keyby = .(patient_id, end.fu)]

dt.tv[, postop.readmit.cohort := start.readmit> 0 & tstop <= final.date.readmit & end.readmit <= 90 & any.op.readmit == T]


dt.tv[,event.readmit :=0]
dt.tv[(postop.readmit.cohort) & 
      emergency_readmitdate  == tstop & 
      COVIDpositivedate != tstop &  
      COVIDreadmission == F & 
      readmit.end > study.start, event.readmit := 1] # COVID isn't a competing risk its an exposure for non covid readmission
dt.tv[(postop.readmit.cohort) & date_death_ons == tstop & event.readmit != 1, event.readmit := 2]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

# Final cohort----

procedures.sub <- c('Colectomy','Cholecystectomy',
                    'HipReplacement','KneeReplacement')
covariates <- c(procedures,procedures.sub,'sex','age.cat','bmi.cat','imd5','wave','LOS.bin',
                'vaccination.status.factor','Current.Cancer','Emergency','Charlson','Charl12','recentCOVID','previousCOVID','postcovid','region',
                'emergency_readmit_primary_diagnosis','primary_diagnosis','death_underlying_cause_ons','admission_method','days_in_critical_care')

drop.vars <- names(dt.tv)[!(names(dt.tv) %in% c(covariates, 'patient_id', 'tstart','tstop','start','end','event','study.start','postop.covid.cohort','end.fu','died','gp.end','discharge.date','age',
                                                'postop.readmit.cohort','postcovid.VTE.cohort','postop.los.cohort','event.VTE','event.readmit','date_death_ons',
                                                'los.end','start.readmit','end.readmit','COVIDpositivedate','COVIDreadmission','emergency_readmitdate','post.VTE.date',
                                                'any.op','any.op.COVID','any.op.readmit','any.op.VTE','admit.date','discharged','final.date','final.date.readmit','final.date.VTE'))]

dt.tv[,(drop.vars) := NULL]


dt.tv <- dt.tv[any.op == T & start >= 0 & tstop <= end.fu,] # Need to start follow up on day after operation as can't identify order when events on same day
arrow::write_feather(dt.tv, sink = here::here("output","cohort_long.feather"))
