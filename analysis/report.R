#detach("package:here", unload = TRUE)
#setwd("C:\\Users\\mczcjc\\Documents\\GitHub\\PostOpCovid")
#library(here)
index_date <- data.table::as.IDate("2020-02-01")
dt <- data.table::fread( here::here("output", "input.csv"))


##########################################################
# Reshape data to long time varying cohort per procedure
#########################################################

### Time splits
fixed <- c('patient_id','dob','sex','bmi', 'region', 'imd','date_death_ons')

time.cols <- c(paste0("covid_vaccine_dates_",1:3),c(names(dt)[grep("^pre",names(dt))]))

proc.tval.stubs <- c('_admission_method','_primary_diagnosis',
'_days_in_critical_care',
'_case_category','_emergency_readmit_primary_diagnosis')

proc.time.stubs <- c('_date_admitted','_date_discharged',
           '_date',
           '_emergency_readmit_date_admitted',
           '_VTE_HES_date_admitted',
           '_VTE_GP_date','_anticoagulation_prescriptions_date')

procedures <- c('LeftHemicolectomy','RightHemicolectomy','TotalColectomy','RectalResection')

proc.time.cols <- paste(rep(procedures,each = length(proc.time.stubs)),proc.time.stubs, sep ="")
proc.tval.cols <- paste(rep(procedures,each = length(proc.tval.stubs)),proc.tval.stubs, sep ="")

dt[,(proc.time.cols) := lapply(.SD,as.numeric), .SDcols = proc.time.cols]
dt[,(time.cols) := lapply(.SD,as.numeric), .SDcols = time.cols]
dt[,date_death_ons := as.numeric(date_death_ons)]

dt[,no.op := rowSums(!is.na(.SD)),.SDcols = c(paste(procedures,"_date_admitted",sep = ""))]
dt[,max.date := lapply(.SD,max, na.rm = T),
   by = patient_id,
   .SDcols = c('date_death_ons',paste(procedures,"_date_discharged",sep = ""))]

dt[is.na(date_death_ons), max.date := max.date + 90]
dt <- dt[no.op > 0,]

dt.fixed <- dt[,.SD, .SDcols = c(fixed,'max.date')]

dt.tv <- survival::tmerge(dt.fixed,dt.fixed,id = patient_id, end.fu = event(max.date) )

dt.times <- dt[,.SD, .SDcols = c('patient_id',time.cols, proc.time.cols)]
dt.tv.values <- dt[,.SD, .SDcols = c('patient_id',paste(procedures,"_date_admitted",sep = ""), proc.tval.cols)]

for (i in c(time.cols,proc.time.cols)) {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.times,",
                           "id = patient_id,", 
                           i," = event(",i,")))")))
}

data.table::setkey(data.table::setDT(dt.tv),patient_id, tstart, tstop)

dt.tv <- dt.tv.values[is.finite(LeftHemicolectomy_date_admitted),.SD,
                      .SDcols =  c('patient_id','LeftHemicolectomy_date_admitted',
                                   paste0("LeftHemicolectomy",proc.tval.stubs))][dt.tv,
                                                                                 on = c(patient_id = "patient_id",LeftHemicolectomy_date_admitted = "tstart")]
data.table::setnames(dt.tv,"LeftHemicolectomy_date_admitted","tstart")
dt.tv[i.LeftHemicolectomy_date_admitted == 1, LeftHemicolectomy_date_admitted := tstart]

dt.tv <- dt.tv.values[is.finite(RightHemicolectomy_date_admitted),.SD,
                      .SDcols = c('patient_id','RightHemicolectomy_date_admitted',
                                  paste0("RightHemicolectomy",proc.tval.stubs))][dt.tv,
                                                                                 on = c(patient_id = "patient_id",RightHemicolectomy_date_admitted = "tstart")]
data.table::setnames(dt.tv,"RightHemicolectomy_date_admitted","tstart")
dt.tv[i.RightHemicolectomy_date_admitted == 1, RightHemicolectomy_date_admitted :=tstart]

dt.tv <- dt.tv.values[is.finite(TotalColectomy_date_admitted),.SD,
                      .SDcols = c('patient_id','TotalColectomy_date_admitted', 
                                  paste0("TotalColectomy",proc.tval.stubs))][dt.tv,
                                                                             on = c(patient_id = "patient_id",TotalColectomy_date_admitted = "tstart")]
data.table::setnames(dt.tv,"TotalColectomy_date_admitted","tstart")
dt.tv[i.TotalColectomy_date_admitted==1,TotalColectomy_date_admitted := tstart]

dt.tv <- dt.tv.values[is.finite(RectalResection_date_admitted),.SD,
                      .SDcols = c('patient_id','RectalResection_date_admitted',
                                  paste0("RectalResection",proc.tval.stubs))][dt.tv,
                                                                              on = c(patient_id = "patient_id",RectalResection_date_admitted = "tstart")]
data.table::setnames(dt.tv,"RectalResection_date_admitted","tstart")
dt.tv[i.RectalResection_date_admitted==1, RectalResection_date_admitted:= tstart]

dt.tv[,c('i.LeftHemicolectomy_date_admitted','i.RightHemicolectomy_date_admitted','i.RectalResection_date_admitted','i.TotalColectomy_date_admitted'):=NULL]
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

admission.dates <- c('admit.date','discharge.date')
dt.tv[,`:=`(admit.date = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_date_admitted")]

dt.tv[LeftHemicolectomy_date_discharged==1 | RightHemicolectomy_date_discharged==1 | RectalResection_date_discharged==1 | TotalColectomy_date_discharged==1, discharge.date := tstop]
dt.tv[, (admission.dates) := lapply(.SD, data.table::nafill, type = "locf"), by = patient_id, .SDcols = admission.dates]
dt.tv[admit.date > discharge.date, (admission.dates) := NA]


dt.tv[,op.type := max(data.table::fifelse(admit.date == LeftHemicolectomy_date_admitted, 'LeftHemicolectomy',
                                      data.table::fifelse(admit.date == RightHemicolectomy_date_admitted, 'RightHemicolectomy',  
                                                          data.table::fifelse(admit.date == TotalColectomy_date_admitted, 'TotalColectomy',
                                                                              data.table::fifelse(admit.date == RectalResection_date_admitted, 'RectalResection','')
                                                          ))),na.rm = T), by = .(patient_id, admit.date)]


data.table::setkey(dt.tv,"patient_id","tstart","tstop")


dt.tv[,`:=`(admission_method = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_admission_method")]
dt.tv[,`:=`(primary_diagnosis = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_primary_diagnosis")]
dt.tv[,`:=`(days_in_critical_care = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_days_in_critical_care")]
dt.tv[,`:=`(case_category = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_case_category")]
dt.tv[,`:=`(emergency_readmit_primary_diagnosis = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_emergency_readmit_primary_diagnosis")]

dt.tv[,(proc.tval.cols) := NULL]

#These flags from tmerge mark event at end of episode
dt.tv[,`:=`(discharged = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_date_discharged")]
dt.tv[,`:=`(COVIDpositive = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_date")]
dt.tv[,`:=`(emergency_readmit = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_emergency_readmit_date_admitted")]
dt.tv[,`:=`(VTE_HES = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_VTE_HES_date_admitted")]
dt.tv[,`:=`(VTE_GP = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_VTE_GP_date")]
dt.tv[,`:=`(anticoagulation_prescriptions = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_anticoagulation_prescriptions_date")]

dt.tv[,(proc.time.cols) := NULL]


###############################
# Pre operative risk factors----
##############################
dt.tv[,age := floor((tstart - as.numeric(as.Date(paste0(dob,'-15'))))/365.25)]
dt.tv[,age.cat := cut(age, breaks = c(18,50,70,80,90),ordered_result = T , right = T, include.lowest = T)]


### Calculate Charlson index at time of operation
comorb.cols <- c(names(dt)[grep("^pre",names(dt))])

data.table::setkey(dt.tv,"patient_id","tstart","tstop")

## tmerge flags end of episode prior to diagnosis, so need to count from subsequent episodes
dt.tv[,(comorb.cols) := lapply(.SD, function(x) cumsum(x) - x), by = patient_id, .SDcols = comorb.cols] 

dt.tv[, Charlson := (pre_MI_GP) +
     (pre_CCF_GP) +
     (pre_PVD_GP) +
     (pre_Stroke_GP) +
     (pre_Dementia_GP) +
     (pre_Respiratory_GP) +
     (pre_RA_SLE_Psoriasis_GP) +
     (pre_Ulcer_or_bleed_GP) +
     (pre_all_liver_GP) + 
     (pre_Cirrhosis_GP)*2 + # counted in all_liver_GP too
     (pre_all_diabetes_GP) +
     (pre_Diabetic_Complications_GP) + # counted in diabetes too
     (pre_Other_Neurology_GP)*2 +
     ((pre_CKD_3_5_GP) | (pre_Renal_GP))*2 +
     (pre_Non_Haematology_malignancy_GP)*2 +
     (pre_Haematology_malignancy_GP)*2 +
     (pre_Metastases_GP)*6 +
     (pre_HIV_GP)*6]  

## Operation type


## Define cancer operations
#dt[,Cancer.Surgery := !is.na(surgery_cancer)] ## TODO  make this more specific to cancer related to operation type
#dt[is.na(Cancer.Surgery), Cancer.Surgery := F]

### Define elective or emergency operations
dt.tv[,Emergency := substr(as.character(admission_method),1,1)=="2"]
dt.tv[is.na(Emergency), Emergency := F]

## Define vaccination status - 14 days post date as effective



##############################
#Post operative outcomes----
##############################

### Length of stay----


### Post operative VTE----
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

VTE.cols <- c('VTE_GP', 'VTE_HES', 'anticoagulation_prescriptions')
dt.tv[,(VTE.cols) := lapply(.SD, function(x) cumsum(ifelse(is.na(x),0,x))), by = patient_id, .SDcols = VTE.cols] 

dt.tv[, post.VTE := (VTE_GP==1 | VTE_HES==1) & 
        anticoagulation_prescriptions == 1 & tstop <= end.fu]  # events flagged at end of episode
dt.tv[post.VTE == T, post.VTE.date := tstop]
dt.tv[,post.VTE.date := min(post.VTE.date), by = patient_id]


### Post operative Covid-19----

### Readmissions----

## Define types of emergency readmissions


### Mortality----
dt.tv[,died := tstop == date_death_ons]
## Cause pf death


#######################
#Crude survival plots----
######################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

dt.tv[end.fu ==1, end90 := tstop]
dt.tv[, end90 := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = 'end90']
dt.tv[, start.spell := min(admit.date), by = .(patient_id, end90) ]

dt.tv[, `:=`(start = tstart - start.spell,
             end = tstop - start.spell)]

crude.surv <- survival::survfit(survival::Surv(start,end,died) ~ op.type, data = dt.tv, id = patient_id)
plot_surv <- survminer::ggsurvplot(crude.surv, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_surv$plot,filename = "plot_surv.png", path=here::here("output"))


crude.los <- survival::survfit(survival::Surv(tstart,tstop,discharged) ~ op.type, data = dt.tv[is.finite(admit.date)], id = patient_id)
plot_los <-survminer::ggsurvplot(crude.los, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_los$plot,filename = "plot_los.png", path=here::here("output"))



##############################
# Post operative COVID risk
#############################

##TODO need to see model output to determine appropriate model fit.

post.op.covid.model <- survival::coxph(survival::Surv(tstart,tstop,COVID_positive) ~ op_type + age.cat + sex + bmi + region +  op_type + Emergency + Cancer.Surgery, id = patient_id, data = dt.tv)
print(xtable::xtable(post.op.covid.model), type = 'html', file = here::here("output","post_op_covid_model.html"))
