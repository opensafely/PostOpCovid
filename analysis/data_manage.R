#detach("package:here", unload = TRUE)
#setwd("C:\\Users\\mczcjc\\Documents\\GitHub\\PostOpCovid")
#library(here)
#detach("package:here", unload = TRUE)
#setwd("P:\\GitHub\\PostOpCovid")
#library(here)
library(data.table)
index_date <- data.table::as.IDate("2020-02-01")
dt <- data.table::fread(here::here("output", "input.csv"))
#########################
# Basic counts and descriptions----
#############################
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

dt[,dateofbirth := (data.table::as.IDate(paste0(dob,'-15')))]
dt[dereg_date != "",gp.end := data.table::as.IDate(paste0(dereg_date,'-15'))]
dt[, imd := as.numeric(imd)]
dt[, imd5 := cut(imd, breaks = seq(0,33000,33000/5),  include.lowest = T, ordered_result = F)]

####################################################################

summary(dt)

lapply(paste0(procedures,"_date_admitted"), function(x) dt[is.finite(get(x)),.N])

dt[,(paste("admit.wave.",procedures, sep ="")) := lapply(.SD, function(x) cut(as.numeric(x), breaks = c(as.numeric(data.table::as.IDate("2020-01-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2020-09-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-05-01", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2021-12-31", format = "%Y-%m-%d")),
                                                                                                        as.numeric(data.table::as.IDate("2022-05-01", format = "%Y-%m-%d"))),
                                                                              labels = c("Wave_1","Wave_2","Wave_3","Wave_4"),
                                                                              ordered = T)), 
                                                         .SDcols = c(paste(procedures,"_date_admitted", sep =""))]
for(x in procedures) {
dt[, (paste0(x,"post.VTE")) := ((!is.na(.SD[,3]) &  
                        .SD[,3] <= .SD[,1] + 90 & .SD[,3] >= .SD[,1]) | 
                       ((!is.na(.SD[,4]) &  
                           .SD[,4] <= .SD[,1] + 90 & .SD[,4] >= .SD[,1]))) & 
        (!is.na(.SD[,5]) &  
           .SD[,5] <= .SD[,1] + 90 & .SD[,5] >= .SD[,1]), 
   .SDcols = paste0(x,c("_date_admitted","_date_discharged","_VTE_GP_date","_VTE_HES_date_admitted","_anticoagulation_prescriptions_date"))]  # events flagged at end of episode
}


demo.waves.tab <- lapply(procedures, function(proc) { 
  t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  {
    cbind(dt[is.finite(get(paste0(proc,'_date_admitted'))) & get(paste0('admit.wave.',proc)) == x,.("Procedures" = .N,
                                                                                                                         "Patients" = length(unique(patient_id)),
                                                                                                                         "Male" = round(mean(sex=='M'),digits = 2),
                                                                                                                         "Age (IQR)" = paste(round(quantile(as.numeric(get(paste0(proc,'_date_admitted')) - as.numeric(dateofbirth))/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                         "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                         "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                         "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_1 >= 14),digits = 2),
                                                                                                                         "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_2 >= 14),digits = 2),
                                                                                                                         "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & get(paste0(proc,'_date_admitted'))  - covid_vaccine_dates_3 >= 14),digits = 2),
                                                                                                                         "Current Cancer"  = round(mean(substr(get(paste0(proc,'_primary_diagnosis')),1,1) =='C'), digits = 2),
                                                                                                                         "Emergency" = round(mean(substr(get(paste0(proc,'_admission_method')),1,1) == "2"),digits = 2),
                                                                                                                         "Length of stay (IQR)" =  paste(round(quantile((get(paste0(proc,'_date_discharged')) - get(paste0(proc,'_date_admitted'))),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                         "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - get(paste0(proc,'_date_admitted')) <= 90),digits = 2),
                                                                                                                         "90 day COVID-19" = round(mean(is.finite(get(paste0(proc,'_date'))) & get(paste0(proc,'_date')) - get(paste0(proc,'_date_admitted')) <= 90  & get(paste0(proc,'_date')) - get(paste0(proc,'_date_admitted')) >=0),digits = 2),
                                                                                                                         "90 day VTE" = round(mean(get(paste0(proc,'post.VTE')), na.rm = T),digits = 2))],
                               t(dt[is.finite(get(paste0(proc,'_date_admitted'))) & get(paste0('admit.wave.',proc)) == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))]),
                               t(dt[get(paste0('admit.wave.',proc)) == x,.N, by = .(.grp = eval(parse(text = paste0(proc,'_primary_diagnosis'))))][order(-N), do.call(paste,c(.SD, sep = ": "))][1:10])                              
                               )}))))})

demo.waves.tab

lapply(1:length(procedures), function(i) print(xtable::xtable(demo.waves.tab[[i]]), type = 'html', here::here("output",paste0("table1",procedures[i],".html"))))

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
proc.time.stubs.start <- c('_date_admitted','_Trauma_HES_date_admitted', 
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
dt[,max.date := lapply(.SD,max, na.rm = T),
   by = patient_id,
   .SDcols = c(paste0(procedures,"_end_fu"))]

dt[is.finite(gp.end) & max.date > gp.end, max.date := gp.end]
dt[!is.finite(max.date), max.date := as.numeric(data.table::as.IDate('2022-02-01'))]

# Data for long cohort table with variables that are fixed at baseline
dt.fixed <- dt[,.SD, .SDcols = c(fixed,'max.date')] # 
dt.tv <- survival::tmerge(dt.fixed,dt.fixed,id = patient_id, end = event(max.date) ) # set survival dataset with final follow up date per patient

# Data for long cohort table with variables defining events
dt.times <- dt[,.SD, .SDcols = c('patient_id',
                                 time.cols,
                                 "gp.end",
                                 proc.time.cols.start, 
                                 proc.time.cols.end,
                                 paste(procedures,"_end_fu",sep = ""))]
dt.times[, gp.end := as.numeric(gp.end)]

# Data for long cohort table with variables that are time varying
dt.tv.values <- dt[,.SD, .SDcols = c('patient_id',
                                     paste(procedures,"_date_admitted",sep = ""), ## Still needed here to define admission date for the values in merging data
                                     paste(procedures,"_emergency_readmit_date_admitted",sep = ""), 
                                     proc.tval.cols)]

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

# Set non event times to missing rather than zero (really only maters for tmerge event)
data.table::setDT(dt.tv)
for (i in c(proc.time.cols.start,
            proc.time.cols.end,
            time.cols,
            proc.tval.cols,
            paste0(procedures,"_end_fu")))  {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = dt.tv[",i,"==0, ",i,":=NA])")))
}

#######
# Copy times across all records within patient and procedure----
#######

# Keep resorting to avoid incorrect copying
data.table::setkey(dt.tv,patient_id, tstart, tstop)


# Copy gp end of follow up across all patient time
dt.tv[, gp.end := max(gp.end, na.rm = T), keyby = patient_id]
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
dt.tv[, (c('discharge.date','end.fu')) := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = c('discharge.date','end.fu')]
dt.tv[discharge.date > end.fu, discharge.date := NA]
dt.tv[tstart >= discharge.date, discharge.date := NA]

## Start of follow up (each enter study) per patient and procedure
data.table::setkey(dt.tv,patient_id, tstart, tstop)
dt.tv[,study.start := min(admit.date, na.rm = T), by = .(patient_id,end.fu)]
dt.tv[!is.finite(study.start), study.start := NA]

## Roll start of procedure periods forward to define beginning of each post procedure period. Include discharge dates and end of fu dates to work out data to drop later
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, (admission.dates) := lapply(.SD, data.table::nafill, type = "locf"), by = patient_id, .SDcols = admission.dates]

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
dt.tv[!is.finite(discharge.date), discharge.date := end.fu]
dt.tv[!is.finite(study.start), study.start := NA]

##########
## Defining operation for each post op period--------------
######
## Assumption can't have two colonic resections on same day
for (i in 1:length(procedures)) dt.tv[study.start == get(paste0(procedures[i],"_date_admitted")),op.type := procedures[i]]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
id_change = dt.tv[, c(TRUE, patient_id[-1] != patient_id[-.N] | end.fu[-1] != end.fu[-.N])]
dt.tv[, op.type := op.type[cummax(((!is.na(x)) | id_change) * .I)]]

dt.tv[is.infinite(op.type), op.type := NA]

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
for(stub in proc.tval.stubs) {dt.tv[,(gsub("*_HES_binary_flag$","",gsub("^_*","",proc.tval.stubs))) :=
                                      data.table::fcoalesce(.SD), .SDcols = paste0(procedures,stub)]}
dt.tv[,Current.Cancer := substr(primary_diagnosis,1,1) =='C']
dt.tv[is.na(Current.Cancer), Current.Cancer := F]
dt.tv[,(proc.tval.cols) := NULL]

## Coalescing outcome event variables into continuous record
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,COVIDpositivedate  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(COVIDpositivedate), COVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,recentCOVIDpositivedate  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_recent_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(recentCOVIDpositivedate), recentCOVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,previousCOVIDpositivedate  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_previous_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(previousCOVIDpositivedate), previousCOVIDpositivedate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,emergency_readmitdate  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_emergency_readmit_date_admitted"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(emergency_readmitdate), emergency_readmitdate := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,VTE_HES  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_VTE_HES_date_admitted"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(VTE_HES), VTE_HES := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,VTE_GP  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_VTE_GP_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(VTE_GP), VTE_GP := NA]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,anticoagulation_prescriptions  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_anticoagulation_prescriptions_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(anticoagulation_prescriptions), anticoagulation_prescriptions := NA]

#data.table::setkey(dt.tv,patient_id,tstart,tstop)
#dt.tv[,Trauma  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_Trauma_HES_date_admitted"), by = .(patient_id, end.fu)]
#dt.tv[!is.finite(Trauma), Trauma := NA]

dt.tv[,(c(proc.time.cols.start,proc.time.cols.end)) := NULL]


###############################
# Pre operative risk factors----
##############################
dt.tv[,age := floor((tstart - as.numeric(as.Date(paste0(dob,'-15'))))/365.25)]
dt.tv[,age.cat := cut(age, breaks = c(1,50,70,80,90),ordered_result = T , right = T, include.lowest = T)]
dt.tv[, imd5 := cut(imd, breaks = seq(0,33000,33000/5),  include.lowest = T, ordered_result = F)]
dt.tv[, bmi.cat := cut(bmi, breaks = c(0,18,24,29,100),  include.lowest = T, ordered_result = F)]


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

dt.tv[,Charl12 := cut(Charlson, breaks = c(0,1,2,100), labels = c("None","Single","Multiple or Severe"), ordered = F)]

## Operation type

#dt.tv[Trauma == F & op.type == 'FractureProcedure',op.type := NA]

## Define cancer operations
#dt[,Cancer.Surgery := !is.na(surgery_cancer)] ## TODO  make this more specific to cancer related to operation type
#dt[is.na(Cancer.Surgery), Cancer.Surgery := F]

### Define elective or emergency operations
dt.tv[,Emergency := substr(as.character(admission_method),1,1)=="2"]
dt.tv[is.na(Emergency), Emergency := F]

## Define vaccination status - 14 days post date as effective
dt.tv[, vaccination.status := is.finite(covid_vaccine_dates_1) + is.finite(covid_vaccine_dates_2) + is.finite(covid_vaccine_dates_3)]
dt.tv[,vaccination.status.factor := factor(vaccination.status, ordered = F)]

##############################
#Post operative outcomes----
##############################

### Length of stay----
dt.tv[,discharged := is.finite(discharge.date) & discharge.date == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,los.end := min(as.numeric(discharge.date), na.rm = T), by = .(patient_id, end.fu) ]

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
dt.tv[,post.VTE.date := min(as.numeric(post.VTE.date), na.rm = T), by = patient_id]
dt.tv[,postVTEany := cumsum(post.VTE), by = .(patient_id, end.fu)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,VTE.end := min(as.numeric(post.VTE.date), na.rm = T), by = .(patient_id, end.fu) ]


### Post operative Covid-19----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,COVIDpositive := is.finite(COVIDpositivedate) & COVIDpositivedate == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,postcovid := cumsum(COVIDpositive), by = .(patient_id, end.fu)]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,covid.end := min(as.numeric(COVIDpositivedate), na.rm = T), by = .(patient_id, end.fu) ]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,recentCOVID := is.finite(recentCOVIDpositivedate) ]
dt.tv[,recentCOVID := max(as.numeric(recentCOVID), na.rm = T), keyby = .(patient_id, end.fu) ]
dt.tv[!is.finite(recentCOVID), recentCOVID := 0]
dt.tv[,previousCOVID := is.finite(previousCOVIDpositivedate) ]
dt.tv[,previousCOVID := max(as.numeric(previousCOVID), na.rm = T), keyby = .(patient_id, end.fu) ]
dt.tv[!is.finite(previousCOVID), previousCOVID := 0]

### Readmissions----
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,emergency_readmit  := is.finite(emergency_readmitdate) & emergency_readmitdate == tstop]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[,readmit.end := min(as.numeric(emergency_readmitdate), na.rm = T), by = .(patient_id, end.fu) ]

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
dt.tv[start >= 0,discharge.start := min(discharge.date, na.rm = T), by = .(patient_id, end.fu)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)

save(dt.tv, file = here::here("output","cohort_long.RData"))