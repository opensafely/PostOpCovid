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
#########################
# Basic counts and descriptions----
#############################
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
procs <- paste0(rep(procedures,each = 5),"_",1:5)

dt[,dateofbirth := (data.table::as.IDate(paste0(dob,'-15')))]
dt[dereg_date != "",gp.end := data.table::as.IDate(paste0(dereg_date,'-15'))]
dt[, imd := as.numeric(imd)]
dt[, imd5 := cut(imd, breaks = seq(-1,33000,33000/5),  include.lowest = T, ordered_result = F)]

####################################################################
# Multiple operations per row. Reshape to long format
####################################################################

#Clean invalid admission sequences
op.admit.vars <- sort(paste0(outer(procedures,paste0("_",1:5),paste0),'_date_admitted'))
op.discharge.vars <- sort(paste0(outer(procedures,paste0("_",1:5),paste0),'_date_discharged'))
for (i in 1:length(op.admit.vars)) {
  dt[get(op.admit.vars[i])>get(op.discharge.vars[i]), 
     (names(dt)[grepl(pattern = paste0(gsub(x =op.admit.vars[i],
                                            pattern = '_date_admitted',
                                            replacement = ""),'*'),
                      x = names(dt))]) := NA]
}
####### Reshape repeated procedures within a specialty
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

#summary(dt)
dt[, keep := F]
dt[!is.na(dereg_date), dereg_date := data.table::as.IDate(paste0(dereg_date,"-30"), format = "%Y-%m-%d")]
for (x in paste0(procedures,"_date_admitted")) { dt[is.na(dereg_date) | x <= dereg_date, keep := T] }
dt <- dt[keep == T,]
lapply(paste0(procedures,"_date_admitted"), function(x) dt[is.finite(get(x)),.N])
# ? reshape each procedure long to get counts, but not unique per patient then

dt[,(paste("admit.wave.",procedures, sep ="")) := lapply(.SD, function(x) cut(as.numeric(x), breaks = c(as.numeric(data.table::as.IDate("2020-01-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2020-09-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-05-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-12-31", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2022-05-01", format = "%Y-%m-%d"))),
                                                                              labels = c("Wave_1","Wave_2","Wave_3","Wave_4"),
                                                                              ordered = T)), 
                                                         .SDcols = c(paste0(procedures,"_date_admitted"))]
for(x in procedures) {
dt[, (paste0(x,"post.VTE")) := ((!is.na(.SD[,3]) &  
                        .SD[,3] <= .SD[,1] + 90 & .SD[,3] >= .SD[,1]) | 
                       ((!is.na(.SD[,4]) &  
                           .SD[,4] <= .SD[,1] + 90 & .SD[,4] >= .SD[,1]))) & 
        (!is.na(.SD[,5]) &  
           .SD[,5] <= .SD[,1] + 90 & .SD[,5] >= .SD[,1]), 
   .SDcols = paste0(procedures,c("_date_admitted","_date_discharged","_VTE_GP_date","_VTE_HES_date_admitted","_anticoagulation_prescriptions_date"))]  # events flagged at end of episode
}

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
##########################################################
# Reshape data to long time varying cohort per procedure ----
#########################################################

### Time splits
fixed <- c('patient_id','dob','sex','bmi' ,'region', 'imd','date_death_ons')

time.cols <- c(paste0("covid_vaccine_dates_",1:3),c(names(dt)[grep("^pre",names(dt))])) 

##Admission exposure information needs to start from beginning of row time period


proc.tval.stubs <- c('_admission_method',
                     '_primary_diagnosis',
                     '_days_in_critical_care',
                     '_case_category',
                     '_recent_case_category',
                     '_previous_case_category',
                     '_emergency_readmit_primary_diagnosis',
                     '_HipReplacement_HES_binary_flag',
                     '_KneeReplacement_HES_binary_flag',
                     '_Cholecystectomy_HES_binary_flag',
                     '_Colectomy_HES_binary_flag')
proc.tval.cols <- paste(rep(procedures,each = length(proc.tval.stubs)),proc.tval.stubs, sep ="")

##Exposure so need to be flagged at start of row time period
proc.time.stubs.start <- c('_date_admitted', 
                           '_recent_date',
                           '_previous_date')
proc.time.cols.start <- paste(rep(procedures,each = length(proc.time.stubs.start)),proc.time.stubs.start, sep ="")

##Outcomes so need to be flagged at end of row time period
proc.time.stubs.end <- c('_date_discharged',
                         '_emergency_readmit_date_admitted',
                           '_date',
                           '_VTE_HES_date_admitted',
                           '_VTE_GP_date',
                         '_anticoagulation_prescriptions_date')
proc.time.cols.end <- paste(rep(procedures,each = length(proc.time.stubs.end)),proc.time.stubs.end, sep ="")

# Dates as numeric for the tmerge
dt[,(c(proc.time.cols.start,proc.time.cols.end)) := lapply(.SD,as.numeric), .SDcols = c(proc.time.cols.start,proc.time.cols.end)]
dt[,(time.cols) := lapply(.SD,as.numeric), .SDcols = time.cols]
dt[,date_death_ons := as.numeric(date_death_ons)]
dt[,gp.end := as.numeric(gp.end)]

# Variable for earliest of 90 days post procedure or death for end of follow up time
dt[,(paste0(procedures,"_end_fu")) := lapply(.SD, function(x) data.table::fifelse(is.finite(date_death_ons) & x+90 > date_death_ons, date_death_ons,x+90)),
   .SDcols = paste0(procedures, '_date_discharged')]
dt[,(paste0(procedures,"_end_fu")) := lapply(.SD, function(x) data.table::fifelse(is.finite(gp.end) & x > gp.end, gp.end,x)),
   .SDcols = paste0(procedures, '_end_fu')]

# Check of how many operations per patients
dt[,no.op := rowSums(!is.na(.SD)),.SDcols = c(paste(procedures,"_date_admitted",sep = ""))]
dt[,summary(no.op)]

# For tmerge define final date for per patient taking into account final patient in GP database and currrent study end date

max.grp.col_(dt = 'dt', max.var.name = 'max.date', aggregate.cols = paste0(procedures,"_end_fu"), id.vars = 'patient_id') 


dt[is.finite(gp.end) & max.date > gp.end, max.date := gp.end]
dt[!is.finite(max.date), max.date := as.numeric(data.table::as.IDate('2022-02-01'))]


# Data for long cohort table with variables that are fixed at baseline
dt.fixed <- unique(dt[,.SD, .SDcols = c(fixed,'max.date')]) # 
dt.tv <- survival::tmerge(dt.fixed,dt.fixed,id = patient_id, end = event(max.date) ) # set survival dataset with final follow up date per patient
rm(dt.fixed)
# Data for long cohort table with variables defining events
dt.times <- dt[,.SD, .SDcols = c('patient_id',
                                 time.cols,
                                 'gp.end',
                                 proc.time.cols.start, 
                                 proc.time.cols.end,
                                 paste(procedures,"_end_fu",sep = ""))]
dt.times[, gp.end := as.numeric(gp.end)]

# Data for long cohort table with variables that are time varying
dt.tv.values <- dt[,.SD, .SDcols = c('patient_id',
                                     paste(procedures,"_date_admitted",sep = ""), ## Still needed here to define admission date for the values in merging data
                                     paste(procedures,"_emergency_readmit_date_admitted",sep = ""), 
                                     proc.tval.cols)]
rm(dt)
# tmerge events defining end of row outcome events
for (i in c(proc.time.cols.end,
            paste0(procedures,"_end_fu"),
            "gp.end"))  {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.times,",
                           "id = patient_id,",
                           i," = event(",i,",",i,")))")))
}

# tmerge events defining start of row exposure events
for (i in c(proc.time.cols.start,
            time.cols)) {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.times,",
                           "id = patient_id,",
                           i," = tdc(",i,",",i,",NA)))")))
}

# tmerge values defining start of row exposure values
for (proc in procedures) {
  for (val in proc.tval.stubs) {
    time.var = '_date_admitted'
    if (val == '_emergency_readmit_primary_diagnosis') time.var = '_emergency_readmit_date_admitted'
       eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.tv.values,",
                             "id = patient_id,", 
                             paste0(proc,val)," = tdc(",paste0(proc,time.var),",",paste0(proc,val),",NA)))")))
  }
}
rm(dt.times)

# Set non event times to missing rather than zero (really only matters for tmerge event)
data.table::setDT(dt.tv)
for (i in c(proc.time.cols.end,
            paste0(procedures,"_end_fu"),
            "gp.end"))  {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = dt.tv[",i,"==0, ",i,":=NA])")))
}

#######
# Copy times across all records within patient and procedure----
#######

# Keep resorting to avoid incorrect copying
data.table::setkey(dt.tv,patient_id, tstart, tstop)

# Copy gp end of follow up across all patient time

max.grp.col_(dt = 'dt.tv', max.var.name = 'gp.end', aggregate.cols = 'gp.end', id.vars = 'patient_id') 

dt.tv[gp.end == 0, gp.end := Inf]

data.table::setkey(dt.tv,patient_id, tstart, tstop)

admission.dates <- c('admit.date','discharge.date','end.fu') # key dates to define

## Coalesce admission from different procedures to create continuous record, only present when valid
dt.tv[,admit.date := do.call(pmax, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_admitted")]
dt.tv[!is.finite(admit.date), admit.date := NA]
dt.tv[admit.date > tstart, admit.date := NA]

## Coalesce end of follow up from different procedures to create continuous record, only present when valid
dt.tv[,end.fu := do.call(pmin, c(.SD, na.rm = T)), .SDcols = c(paste0(procedures,"_end_fu"))] ## gp.end in end_fu already
dt.tv[!is.finite(end.fu), end.fu := NA]
dt.tv[end.fu <= tstart, end.fu := NA]

## Coalesce discharged date from different procedures to create continuous record, only present when valid
dt.tv[,discharge.date := do.call(pmax, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_discharged")]
dt.tv[!is.finite(discharge.date), discharge.date := NA]
data.table::setkey(dt.tv,patient_id,tstart,tstop)

## Roll end of period dates backwards to end of previous episodes to define each post procedure period
data.table::setkey(dt.tv,patient_id, tstart, tstop)
nocb.roll_(dt = 'dt.tv', ID = 'patient_id', start.DTTM = 'tstart', group = 'c("patient_id")', var.cols = 'c("discharge.date","end.fu")')

dt.tv[, (c('discharge.date','end.fu')) := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = c('discharge.date','end.fu')]
dt.tv[discharge.date > end.fu, discharge.date := NA]
dt.tv[tstart >= discharge.date, discharge.date := NA]

## Start of follow up (each enter study) per patient and procedure
data.table::setkey(dt.tv,patient_id, tstart, tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'study.start',aggregate.cols = 'admit.date',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(study.start), study.start := NA]

## Roll start of procedure periods forward to define beginning of each post procedure period. Include discharge dates and end of fu dates to work out data to drop later
data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv', ID = 'patient_id', start.DTTM = 'tstart', group = 'c("patient_id")', var.cols = 'c(admission.dates)')

## Remove admit & discharge dates from outside admission period. Study start and endfu define total exposure for each observation period.
dt.tv[admit.date > discharge.date | is.na(admit.date) | admit.date > tstart | discharge.date < tstop, c('admit.date','discharge.date') := NA]
dt.tv[is.na(admit.date) , (paste0(procedures,c('_admission_method','_primary_diagnosis',
                                                        '_days_in_critical_care')))  := NA]
dt.tv <- dt.tv[tstop <= end.fu,]

# Pre study start not dropped yet as contains pre operative exposure information

#for (proc in procedures) {
#  cols <-paste0(proc, proc.time.stubs.start,proc.time.stubs.end)
#  eval(parse(text = paste0("dt.tv[",proc,"_date_admitted>",proc,"_date_discharged,(cols) := NA]")))
#}

dt.tv[!is.finite(admit.date), admit.date := NA]
dt.tv[!is.finite(end.fu), end.fu := NA]
dt.tv[!is.finite(discharge.date), discharge.date := NA]
dt.tv[!is.finite(study.start), study.start := NA]

##########
## Defining operation for each post op period--------------
######
## Assumption can't be made that can't have more than procedures on same day

for (i in 1:length(procedures)) dt.tv[, (procedures[i]) := is.finite(get(paste0(procedures[i],"_date_admitted"))) &
                                        is.finite(admit.date) &
                                        admit.date <= get(paste0(procedures[i],"_date_admitted")) & 
                                        (is.finite(end.fu) | 
                                           end.fu >= get(paste0(procedures[i],"_date_admitted")))]
dt.tv[,(procedures) := lapply(.SD, function(x) data.table::fifelse(x==0,NA,x)), .SDcols = c(procedures)]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
locf.roll_(dt = 'dt.tv', ID = 'patient_id', start.DTTM = 'tstart', group = 'c("patient_id","end.fu")', var.cols = 'c(procedures)')
data.table::setkey(dt.tv,patient_id,tstart,tstop)
for (i in 1:length(procedures)) dt.tv[!is.finite(get(procedures[i])), (procedures[i]) := F]

## Coalesce across values for each post op exposure period. Already carried forward in time by tmerge
data.table::setkey(dt.tv,patient_id,tstart,tstop)


proc.tval.stubs <- c('_admission_method',
                     '_primary_diagnosis',
                     '_days_in_critical_care',
                     '_case_category',
                     '_emergency_readmit_primary_diagnosis',
                     '_HipReplacement_HES_binary_flag',
                     '_KneeReplacement_HES_binary_flag',
                     '_Cholecystectomy_HES_binary_flag',
                     '_Colectomy_HES_binary_flag')
for(stub in proc.tval.stubs) {dt.tv[,(gsub("*_HES_binary_flag$","",gsub("^_*","",stub))) :=
                                      data.table::fcoalesce(.SD), .SDcols = paste0(procedures,stub)]}
dt.tv[,Current.Cancer := substr(primary_diagnosis,1,1) =='C']
dt.tv[is.na(Current.Cancer), Current.Cancer := F]
dt.tv[,(proc.tval.cols) := NULL]

sub.procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')
dt.tv[,(sub.procedures) := lapply(.SD, function(x) max(x, na.rm = T)), by = .(patient_id, end.fu), .SDcols = c(sub.procedures)]                

## Coalescing outcome event variables into continuous record
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'COVIDpositivedate',aggregate.cols = paste0(procedures,"_date"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(COVIDpositivedate), COVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'recentCOVIDpositivedate',aggregate.cols = paste0(procedures,"_recent_date"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(recentCOVIDpositivedate), recentCOVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'previousCOVIDpositivedate',aggregate.cols = paste0(procedures,"_previous_date"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(previousCOVIDpositivedate), previousCOVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'emergency_readmitdate',aggregate.cols = paste0(procedures,"_emergency_readmit_date_admitted"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(emergency_readmitdate), emergency_readmitdate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'VTE_HES',aggregate.cols = paste0(procedures,"_VTE_HES_date_admitted"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(VTE_HES), VTE_HES := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'VTE_GP',aggregate.cols = paste0(procedures,"_VTE_GP_date"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(VTE_GP), VTE_GP := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'anticoagulation_prescriptions',aggregate.cols = paste0(procedures,"_anticoagulation_prescriptions_date"),id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(anticoagulation_prescriptions), anticoagulation_prescriptions := NA]

#data.table::setkey(dt.tv,patient_id,tstart,tstop)
#dt.tv[,Trauma  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_Trauma_HES_date_admitted"), by = .(patient_id, end.fu)]
#dt.tv[!is.finite(Trauma), Trauma := NA]

dt.tv[,(c(proc.time.cols.start,proc.time.cols.end)) := NULL]


###############################
# Pre operative risk factors----
##############################
dt.tv[,age := floor((tstart - as.numeric(as.Date(paste0(dob,'-15'))))/365.25)]
max.age <- max(dt.tv$age,na.rm = T)
dt.tv[,age.cat := cut(age, breaks = c(1,50,70,80,max.age),ordered_result = F , right = T, include.lowest = T)]

dt.tv[, imd5 := cut(imd, breaks = seq(0,33000,33000/5),  include.lowest = T, ordered_result = F)]
levels(dt.tv$imd5) <- c(levels(dt.tv$imd5),"Missing")
dt.tv[is.na(imd5) , imd5 := "Missing"]

dt.tv[, bmi.cat := cut(bmi, breaks = c(0,18,24,29,100),  include.lowest = T, ordered_result = F)]
levels(dt.tv$bmi.cat) <- c(levels(dt.tv$bmi.cat),"Missing")
dt.tv[is.na(bmi.cat) , bmi.cat := "Missing"]

### Calculate Charlson index at time of operation - tdc so date present from first recording
comorb.cols <- c(names(dt)[grep("^pre",names(dt))])

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
levels(dt.tv$Charl12) <- c(levels(dt.tv$Charl12),"Missing")
dt.tv[is.na(Charl12) , Charl12 := "Missing"]
## Operation type

#dt.tv[Trauma == F & op.type == 'FractureProcedure',op.type := NA]

## Define cancer operations
#dt[,Cancer.Surgery := !is.na(surgery_cancer)] ## TODO  make this more specific to cancer related to operation type
#dt[is.na(Cancer.Surgery), Cancer.Surgery := F]

### Define elective or emergency operations
dt.tv[,Emergency := substr(as.character(admission_method),1,1)=="2"]
dt.tv[is.na(Emergency), Emergency := F]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'Emergency',aggregate.cols = 'Emergency',id.vars = c("patient_id","end.fu"))

## Define vaccination status - 14 days post date as effective
dt.tv[, vaccination.status := is.finite(covid_vaccine_dates_1) + is.finite(covid_vaccine_dates_2) + is.finite(covid_vaccine_dates_3)]
dt.tv[,vaccination.status.factor := factor(vaccination.status,  ordered = F)]
dt.tv[is.na(vaccination.status.factor), vaccination.status.factor := 0]
##############################
#Post operative outcomes----
##############################

### Length of stay----
dt.tv[,discharged := is.finite(discharge.date) & discharge.date == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'los.end',aggregate.cols = 'discharge.date',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(los.end), los.end := end.fu]

### Post operative VTE----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, anticoagulation_prescriptions := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = 'anticoagulation_prescriptions']

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, post.VTE := ((is.finite(VTE_GP) &  
                        VTE_GP == tstop) | 
                       (is.finite(VTE_HES) & 
                         VTE_HES == tstop)) & 
        is.finite(anticoagulation_prescriptions) &
        tstop <= end.fu]  # events flagged at end of episode
dt.tv[post.VTE == T, post.VTE.date := tstop]
min.grp.col_(dt = 'dt.tv',min.var.name = 'post.VTE.date',aggregate.cols = 'post.VTE.date',id.vars = c("patient_id","end.fu"))
dt.tv[,postVTEany := cumsum(post.VTE), by = .(patient_id, end.fu)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'VTE.end',aggregate.cols = 'post.VTE.date',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(VTE.end), VTE.end := end.fu]

### Post operative Covid-19----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,COVIDpositive := is.finite(COVIDpositivedate) & COVIDpositivedate == tstop]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,postcovid := cumsum(COVIDpositive), by = .(patient_id, end.fu)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, covid.end := COVIDpositivedate]
dt.tv[!is.finite(covid.end), covid.end := end.fu]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,recentCOVID := is.finite(recentCOVIDpositivedate) ]

max.grp.col_(dt = 'dt.tv',max.var.name = 'recentCOVID',aggregate.cols = 'recentCOVID',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(recentCOVID), recentCOVID := 0]

dt.tv[,previousCOVID := is.finite(previousCOVIDpositivedate) ]
max.grp.col_(dt = 'dt.tv',max.var.name = 'previousCOVID',aggregate.cols = 'previousCOVID',id.vars = c("patient_id","end.fu"))
dt.tv[!is.finite(previousCOVID), previousCOVID := 0]

### Readmissions----
dt.tv[,emergency_readmit  := is.finite(emergency_readmitdate) & emergency_readmitdate == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
min.grp.col_(dt = 'dt.tv',min.var.name = 'readmit.end',aggregate.cols = 'emergency_readmitdate',id.vars = c("patient_id","end.fu"))

## Define types of emergency readmissions



### Mortality----
## date_death_ons part of definition of end_fu so will be end of final row when in follow up period
dt.tv[,died := is.finite(date_death_ons) & tstop == date_death_ons]
dt.tv[is.na(died), died := 0]
## Cause of death TODO

data.table::setkey(dt.tv,patient_id,tstart,tstop)
## Pre operation exposures now defined so can drop prior to analysis
dt.tv <- dt.tv[!(tstart < study.start | tstop > end.fu) & !is.na(age.cat),]
dt.tv[, year := data.table::year(data.table::as.IDate(admit.date))]

# Redefine wave in long table
dt.tv[,wave := cut(study.start, breaks = c(as.numeric(data.table::as.IDate("2020-01-01")),
                                           as.numeric(data.table::as.IDate("2020-09-01")),
                                           as.numeric(data.table::as.IDate("2021-05-01")),
                                           as.numeric(data.table::as.IDate("2021-12-31")),
                                           as.numeric(data.table::as.IDate("2022-05-01"))),
                   labels = c("Wave_1","Wave_2","Wave_3","Wave_4"),
                   include.lowest = T,
                   right = T,
                   ordered = F)]

# Restart clock with each procedure
dt.tv[, `:=`(start = tstart - study.start,
             end = tstop - study.start)]
min.grp.col_(dt = 'dt.tv[start >= 0,]',min.var.name = 'discharge.start',aggregate.cols = 'discharge.date',id.vars = c("patient_id","end.fu"))

data.table::setkey(dt.tv,patient_id,tstart,tstop)


############## Define cohorts


### post op COVID cohort
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,final.date := covid.end]
dt.tv[is.finite(readmit.end) & readmit.end < final.date & readmit.end > study.start, final.date := readmit.end]
dt.tv[is.finite(end.fu) & end.fu < final.date, final.date := end.fu]
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date',aggregate.cols = 'final.date',id.vars = c("patient_id","end.fu"))


dt.tv[,event :=0]
dt.tv[COVIDpositivedate == tstop, event := 1]
dt.tv[emergency_readmitdate  == tstop & event != 1, event := 2]
dt.tv[date_death_ons == tstop & event != 1, event := 3]

dt.tv[, postop.covid.cohort := start>=0 & tstop <= final.date & end <= 90]

dt.tv[(postop.covid.cohort) & start ==0  & is.finite(admit.date),any.op := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[is.na(any.op), any.op := F]
dt.tv[, any.op := any.op > 0]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, any.op := cummax(any.op), keyby = .(patient_id, end.fu)]

dt.tv[, postop.covid.cohort := start>=0 & tstop <= final.date & end <= 90 & any.op == T]

### post COVID VTE cohort
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,final.date.VTE := VTE.end]
dt.tv[is.finite(readmit.end) & readmit.end < final.date.VTE, final.date.VTE := readmit.end]
dt.tv[is.finite(end.fu) & end.fu < final.date.VTE, final.date.VTE := end.fu]
min.grp.col_(dt = 'dt.tv',min.var.name = 'final.date.VTE',aggregate.cols = 'final.date.VTE',id.vars = c("patient_id","end.fu"))

dt.tv[,event.VTE :=0]
dt.tv[post.VTE.date == tstop, event.VTE := 1]
dt.tv[emergency_readmitdate  == tstop & event.VTE != 1, event.VTE := 2]
dt.tv[date_death_ons == tstop & event.VTE != 1, event.VTE := 3]

dt.tv[, postcovid.VTE.cohort := start>=0 & tstop <= final.date.VTE & end <= 90]
dt.tv[(postcovid.VTE.cohort) & start ==0  & is.finite(admit.date),any.op.VTE := rowSums(.SD,na.rm =T), .SDcols = c(procedures)]
dt.tv[, any.op := cummax(any.op.VTE), keyby = .(patient_id, end.fu)]

dt.tv[, postcovid.VTE.cohort := start>=0 & tstop <= final.date.VTE & end <= 90 & any.op.VTE == T]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv <- dt.tv[any.op == T & start >=0 & end <=90,]
save(dt.tv, file = here::here("output","cohort_long.RData"))