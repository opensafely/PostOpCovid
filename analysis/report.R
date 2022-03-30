#detach("package:here", unload = TRUE)
#setwd("C:\\Users\\mczcjc\\Documents\\GitHub\\PostOpCovid")
#library(here)
index_date <- data.table::as.IDate("2020-02-01")
dt <- data.table::fread( here::here("output", "input.csv"))
data.table::setkey(dt,"patient_id","surgery_date")
###############################
# Pre operative risk factors----
##############################

dt[,age.cat := cut(age, breaks = c(18,50,70,80,90),ordered_result = T , right = T, include.lowest = T)]

### Calculate Charlson index at time of operation

dt[, Charlson := is.finite(MI_GP) +
     is.finite(CCF_GP) +
     is.finite(PVD_GP) +
     is.finite(Stroke_GP) +
     is.finite(Dementia_GP) +
     is.finite(Chronic_Respiratory_GP) +
     is.finite(RA_SLE_Psoriasis_GP) +
     is.finite(Ulcer_or_bleed_GP) +
     is.finite(all_liver_GP) + 
     is.finite(cirrhosis_GP)*2 + # counted in all_liver_GP too
     is.finite(all_diabetes_GP) +
     is.finite(diabetic_complication_GP) + # counted in diabetes too
     is.finite(other_neuro_GP)*2 +
     (is.finite(CKD_3_5_GP) | is.finite(Renal_GP))*2 +
     is.finite(non_haem_cancer_GP)*2 +
     is.finite(haem_cancer_GP)*2 +
     is.finite(metastatic_cancer_GP)*6 +
     is.finite(HIV_GP)*6]  

## Operation type

# dt[,op_type := as.factor(ifelse(right_hemicolectomy == op.date, "RightHemicolectomy",
#                             ifelse(left_hemicolectomy == op.date, "LeftHemicolectomy",
#                                    ifelse(total_colectomy == op.date, "TotalColectomy",
#                                           ifelse(rectal_resection == op.date,"RectalResection",0)))))]
dt[,op_type := as.factor(ifelse(is.finite(right_hemicolectomy), "RightHemicolectomy",
                                ifelse(is.finite(left_hemicolectomy), "LeftHemicolectomy",
                                       ifelse(is.finite(total_colectomy), "TotalColectomy",
                                              ifelse(is.finite(rectal_resection),"RectalResection",0)))))]


## Define cancer operations
dt[,Cancer.Surgery := !is.na(surgery_cancer)] ## TODO  make this more specific to cancer related to operation type
dt[is.na(Cancer.Surgery), Cancer.Surgery := F]

### Define elective or emergency operations
dt[,Emergency := substr(as.character(surgery_admimeth),1,1)=="2"]
dt[is.na(Emergency), Emergency := F]

##############################
#Post operative outcomes----
##############################

### Post operative VTE----

dt[(is.finite(VTE_GP) | is.finite(VTE_HES)) & 
     Anticoagulant_prescription < surgery_discharge_date + 90,
   post.VTE := min(VTE_GP, VTE_HES, na.rm = T)]

### Time splits

dt[, op.date := min( right_hemicolectomy,left_hemicolectomy, total_colectomy,rectal_resection, na.rm = T)]
#dt[, end.date := data.table::as.IDate("2022-02-01")]
dt[, end.date := surgery_discharge_date + 90]

event.times <- c('surgery_discharge_date', 
                'SARS_CoV_2_test_date', 'VTE_GP', 'date_death_ons', 
                'end.date', 
                'first_emergency_readmission_date','post.VTE')

dt.times <- dt[,.SD, .SDcols = c('patient_id','surgery_date',event.times)]
dt.times[,surgery_date := as.numeric(surgery_date)]
dt.times[,(event.times) := lapply(.SD,as.numeric), .SDcols = event.times]


dt.fixed <- dt[,.(patient_id,age,sex,bmi, region, imd, age.cat,
                  surgery_date, op_type,Emergency, Cancer.Surgery, index_ICU_days,
                  index_primary_diagnosis,first_emergency_readmission_diagnosis)]
####? TODO change to loop
dt.tv <- survival::tmerge(dt.fixed,dt.times,
                 id = patient_id, 
                 end.study = event(end.date))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          start.study = event(surgery_date))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          surgery.discharge = event(surgery_discharge_date))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          death = event(date_death_ons))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          COVID_positive = event(SARS_CoV_2_test_date))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          emergency.readmit = event(first_emergency_readmission_date))

dt.tv <- survival::tmerge(dt.tv,dt.times,
                          id = patient_id, 
                          VTE90days = event(post.VTE))

dt.tv <- data.table::data.table(dt.tv,key = c('patient_id','tstart'))
dt.tv[, `:=`(tstart = as.numeric(tstart),
             tstop = as.numeric(tstop))]

dt.tv <- dt.tv[surgery_date <= tstart]



dt.tv[, `:=`(tstart = tstart - as.numeric(surgery_date),
             tstop = tstop - as.numeric(surgery_date))]
### Post operative Covid-19----

### Readmissions----

## Define types of emergency readmissions


### Mortality----

## Cause pf death


#######################
#Crude survival plots----
######################
crude.surv <- survival::survfit(survival::Surv(tstart,tstop,death) ~ op_type, data = dt.tv, id = patient_id)
plot_surv <-survminer::ggsurvplot(crude.surv, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_surv$plot,filename = "plot_surv.png", path=here::here("output"))


crude.los <- survival::survfit(survival::Surv(tstart,tstop,surgery.discharge) ~ op_type, data = dt.tv, id = patient_id)
plot_los <-survminer::ggsurvplot(crude.surv, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_los$plot,filename = "plot_los.png", path=here::here("output"))



##############################
# Post operative COVID risk
#############################

##TODO need to see model output to determine appropriate model fit.

post.op.covid.model <- survival::coxph(survival::Surv(tstart,tstop,COVID_positive) ~ op_type + age.cat + sex + bmi + region +  op_type + Emergency + Cancer.Surgery, id = patient_id, data = dt.tv)
print(xtable::xtable(post.op.covid.model), type = 'html', file = here::here("output","post_op_covid_model.html"))
