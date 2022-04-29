#detach("package:here", unload = TRUE)
#setwd("C:\\Users\\mczcjc\\Documents\\GitHub\\PostOpCovid")
#library(here)
#detach("package:here", unload = TRUE)
#setwd("P:\\GitHub\\PostOpCovid")
#library(here)
library(data.table)
index_date <- data.table::as.IDate("2020-02-01")
dt <- data.table::fread( here::here("output", "input.csv"))
#########################
# Basic counts and descriptions
#############################

dt[,dateofbirth := (data.table::as.IDate(paste0(dob,'-15')))]
dt[dereg_date != "",gp.end := data.table::as.IDate(paste0(dereg_date,'-15'))]

###Create logical pseudo dates
# dt[LeftHemicolectomy_date_discharged < LeftHemicolectomy_date_admitted,LeftHemicolectomy_date_discharged := LeftHemicolectomy_date_admitted + sample.int(90,1)]
# dt[RightHemicolectomy_date_admitted>= LeftHemicolectomy_date_admitted &  RightHemicolectomy_date_admitted <LeftHemicolectomy_date_admitted + 90 , `:=`(RightHemicolectomy_date_admitted = NA, RightHemicolectomy_date_discharged = NA)]
# dt[TotalColectomy_VTE_HES_date_admitted>= RightHemicolectomy_date_admitted & TotalColectomy_date_admitted <RightHemicolectomy_date_admitted + 90, `:=`(TotalColectomy_date_admitted = NA, TotalColectomy_date_discharged = NA)]
# dt[RectalResection_date_admitted>= TotalColectomy_date_admitted & RectalResection_date_admitted <TotalColectomy_date_admitted + 90, `:=`(RectalResection_date_admitted = NA, RectalResection_date_discharged = NA)]
# 
# procedures <- c('LeftHemicolectomy','RightHemicolectomy','TotalColectomy','RectalResection')
# 
# cols <- paste0(procedures,c("_date_discharged"))
# dt[, (cols) := lapply(.SD, function(x) fifelse(x > date_death_ons,date_death_ons, x)), .SDcols = cols]
# dt[, (cols) := lapply(.SD, function(x) fifelse(x > gp.end,gp.end, x)), .SDcols = cols]
# 
# 
# dt[,last.proc := pmax(TotalColectomy_date_admitted,LeftHemicolectomy_date_admitted,RightHemicolectomy_date_admitted,RectalResection_date_admitted, na.rm = T)]
# dt[date_death_ons < last.proc , date_death_ons := last.proc + sample.int(90,1), ]
# 
# 
# cols <- paste0(procedures,c("_date_discharged"))
# dt[, (cols) := lapply(.SD, function(x) fifelse(x > date_death_ons,date_death_ons, x)), .SDcols = cols]
# 
# cols <- paste0(procedures,c("_date_admitted"))
# dt[, (cols) := lapply(.SD, function(x) ifelse(x > date_death_ons,NA, x)), .SDcols = cols]
# dt[, (cols) := lapply(.SD, function(x) ifelse(x > gp.end,NA, x)), .SDcols = cols]
# 
# 
# dt[LeftHemicolectomy_date_admitted >LeftHemicolectomy_date_discharged,  LeftHemicolectomy_date_discharged := pmin(LeftHemicolectomy_date_admitted + sample.int(90,1), date_death_ons)]
# dt[RightHemicolectomy_date_admitted >RightHemicolectomy_date_discharged,  RightHemicolectomy_date_discharged := pmin(RightHemicolectomy_date_admitted + sample.int(90,1), date_death_ons)]
# dt[TotalColectomy_date_admitted >TotalColectomy_date_discharged,  TotalColectomy_date_discharged := pmin(TotalColectomy_date_admitted + sample.int(90,1), date_death_ons)]
# dt[RectalResection_date_admitted >RectalResection_date_discharged,  RectalResection_date_discharged := pmin(RectalResection_date_admitted + sample.int(90,1), date_death_ons)]


####################################################################


summary(dt)

dt[is.finite(LeftHemicolectomy_date_admitted),.N]
dt[is.finite(RightHemicolectomy_date_admitted),.N]
dt[is.finite(TotalColectomy_date_admitted),.N]
dt[is.finite(RectalResection_date_admitted),.N]

procedures <- c('LeftHemicolectomy','RightHemicolectomy','TotalColectomy','RectalResection')


dt[,(paste("admit.wave.",procedures, sep ="")) := lapply(.SD, function(x) factor(data.table::fifelse(x <= as.numeric(data.table::as.IDate("2020-09-01")),"Wave_1",
                                           data.table::fifelse(x <= as.numeric(data.table::as.IDate("2021-05-01")),"Wave_2",
                                                               data.table::fifelse(x <= as.numeric(data.table::as.IDate("2021-12-31")),"Wave_3","Wave_4"))), ordered = F)), .SDcols = c(paste(procedures,"_date_admitted", sep =""))]
for(x in procedures) {
dt[, (paste0(x,"post.VTE")) := ((!is.na(.SD[,3]) &  
                        .SD[,3] <= .SD[,1] + 90 & .SD[,3] >= .SD[,1]) | 
                       ((!is.na(.SD[,4]) &  
                           .SD[,4] <= .SD[,1] + 90 & .SD[,4] >= .SD[,1]))) & 
        (!is.na(.SD[,5]) &  
           .SD[,5] <= .SD[,1] + 90 & .SD[,5] >= .SD[,1]), 
   .SDcols = paste0(x,c("_date_admitted","_date_discharged","_VTE_GP_date","_VTE_HES_date_admitted","_anticoagulation_prescriptions_date"))]  # events flagged at end of episode
}

t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  
data.table::rbindlist(list(dt[is.finite(LeftHemicolectomy_date_admitted) & admit.wave.LeftHemicolectomy == x,.("Procedures" = .N,
                                                                                    "Patients" = length(unique(patient_id)),
                                                                                    "Male" = round(mean(sex=='M'),digits = 2),
                                                                                    "Age (IQR)" = paste(round(quantile((LeftHemicolectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                    "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                    "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                    "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_1 >= 14),digits = 2),
                                                                                    "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_2 >= 14),digits = 2),
                                                                                    "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_3 >= 14),digits = 2),
                                                                                    "Current Cancer"  = round(mean(substr(LeftHemicolectomy_primary_diagnosis,1,1) =='C'), digits = 2),
                                                                                    "Emergency" = round(mean(substr(LeftHemicolectomy_admission_method,1,1) == "2"),digits = 2),
                                                                                    "Length of stay (IQR)" =  paste(round(quantile((LeftHemicolectomy_date_discharged - LeftHemicolectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                    "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - LeftHemicolectomy_date_admitted <= 90),digits = 2),
                                                                                    "90 day COVID-19" = round(mean(is.finite(LeftHemicolectomy_date) & LeftHemicolectomy_date - LeftHemicolectomy_date_admitted <= 90  & LeftHemicolectomy_date - LeftHemicolectomy_date_admitted >=0),digits = 2),
                                                                                    "90 day VTE" = round(mean(LeftHemicolectomypost.VTE, na.rm = T),digits = 2))]))))))

t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  
  data.table::rbindlist(list(dt[is.finite(RightHemicolectomy_date_admitted) & admit.wave.RightHemicolectomy == x,.("Procedures" = .N,
                                                                                                                 "Patients" = length(unique(patient_id)),
                                                                                                                 "Male" = round(mean(sex=='M'),digits = 2),
                                                                                                                 "Age (IQR)" = paste(round(quantile((RightHemicolectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_1 >= 14),digits = 2),
                                                                                                                 "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_2 >= 14),digits = 2),
                                                                                                                 "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_3 >= 14),digits = 2),
                                                                                                                 "Current Cancer"  = round(mean(substr(RightHemicolectomy_primary_diagnosis,1,1) =='C'), digits = 2),
                                                                                                                 "Emergency" = round(mean(substr(RightHemicolectomy_admission_method,1,1) == "2"),digits = 2),
                                                                                                                 "Length of stay (IQR)" =  paste(round(quantile((RightHemicolectomy_date_discharged - RightHemicolectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - RightHemicolectomy_date_admitted <= 90),digits = 2),
                                                                                                                 "90 day COVID-19" = round(mean(is.finite(RightHemicolectomy_date) & RightHemicolectomy_date - RightHemicolectomy_date_admitted <= 90  & RightHemicolectomy_date - RightHemicolectomy_date_admitted >=0),digits = 2),
                                                                                                                 "90 day VTE" = round(mean(RightHemicolectomypost.VTE, na.rm = T),digits = 2))]))))))

t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  
  data.table::rbindlist(list(dt[is.finite(TotalColectomy_date_admitted) & admit.wave.TotalColectomy == x,.("Procedures" = .N,
                                                                                                                 "Patients" = length(unique(patient_id)),
                                                                                                                 "Male" = round(mean(sex=='M'),digits = 2),
                                                                                                                 "Age (IQR)" = paste(round(quantile((TotalColectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & TotalColectomy_date_admitted  - covid_vaccine_dates_1 >= 14),digits = 2),
                                                                                                                 "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & TotalColectomy_date_admitted  - covid_vaccine_dates_2 >= 14),digits = 2),
                                                                                                                 "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & TotalColectomy_date_admitted  - covid_vaccine_dates_3 >= 14),digits = 2),
                                                                                                           "Current Cancer"  = round(mean(substr(TotalColectomy_primary_diagnosis,1,1) =='C'), digits = 2),
                                                                                                           "Emergency" = round(mean(substr(TotalColectomy_admission_method,1,1) == "2"),digits = 2),
                                                                                                                 "Length of stay (IQR)" =  paste(round(quantile((TotalColectomy_date_discharged - TotalColectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - TotalColectomy_date_admitted <= 90),digits = 2),
                                                                                                                 "90 day COVID-19" = round(mean(is.finite(TotalColectomy_date) & TotalColectomy_date - TotalColectomy_date_admitted <= 90  & TotalColectomy_date - TotalColectomy_date_admitted >=0),digits = 2),
                                                                                                                 "90 day VTE" = round(mean(TotalColectomypost.VTE, na.rm = T),digits = 2))]))))))

t(data.table::rbindlist((lapply(paste0("Wave_",1:4), function(x)  
  data.table::rbindlist(list(dt[is.finite(RectalResection_date_admitted) & admit.wave.RectalResection == x,.("Procedures" = .N,
                                                                                                                 "Patients" = length(unique(patient_id)),
                                                                                                                 "Male" = round(mean(sex=='M'),digits = 2),
                                                                                                                 "Age (IQR)" = paste(round(quantile((RectalResection_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "BMI (IQR)" = paste(round(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "IMD (IQR)" = paste(round(quantile(imd,c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "1st Vaccination" = round(mean(is.finite(covid_vaccine_dates_1) & RectalResection_date_admitted  - covid_vaccine_dates_1 >= 14),digits = 2),
                                                                                                                 "2nd Vaccination" = round(mean(is.finite(covid_vaccine_dates_2) & RectalResection_date_admitted  - covid_vaccine_dates_2 >= 14),digits = 2),
                                                                                                                 "3rd Vaccination" = round(mean(is.finite(covid_vaccine_dates_3) & RectalResection_date_admitted  - covid_vaccine_dates_3 >= 14),digits = 2),
                                                                                                             "Current Cancer"  = round(mean(substr(RectalResection_primary_diagnosis,1,1) =='C'), digits = 2),
                                                                                                             "Emergency" = round(mean(substr(RectalResection_admission_method,1,1) == "2"),digits = 2),
                                                                                                                 "Length of stay (IQR)" =  paste(round(quantile((RectalResection_date_discharged - RectalResection_date_admitted),c(0.25,0.5,0.75),na.rm = T),digits = 2),collapse = ","),
                                                                                                                 "90 day mortality" = round(mean(is.finite(date_death_ons) & date_death_ons - RectalResection_date_admitted <= 90),digits = 2),
                                                                                                                 "90 day COVID-19" = round(mean(is.finite(RectalResection_date) & RectalResection_date - RectalResection_date_admitted <= 90  & RectalResection_date - RectalResection_date_admitted >=0),digits = 2),
                                                                                                                 "90 day VTE" = round(mean(RectalResectionpost.VTE, na.rm = T),digits = 2))]))))))


lapply(paste0("Wave_",1:4), function(x)  
print(xtable::xtable(rbind(t(rbind(dt[is.finite(LeftHemicolectomy_date_admitted) & admit.wave.LeftHemicolectomy == x,.("Procedures" = .N,
                                                  "Patients" = length(unique(patient_id)),
                                                  "Male" = mean(sex=='M'),
                                                  "Age (IQR)" = paste(quantile((LeftHemicolectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                 "BMI (IQR)" = paste(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                 "IMD (IQR)" = paste(quantile(imd,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                 "1st Vaccination" = mean(is.finite(covid_vaccine_dates_1) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_1 >= 14),
                                                 "2nd Vaccination" = mean(is.finite(covid_vaccine_dates_2) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_2 >= 14),
                                                 "3rd Vaccination" = mean(is.finite(covid_vaccine_dates_3) & LeftHemicolectomy_date_admitted  - covid_vaccine_dates_3 >= 14),
                                                  "Emergency" = mean(substr(LeftHemicolectomy_admission_method,1,1) == "2"),
                                                  "Length of stay (IQR)" =  paste(quantile((LeftHemicolectomy_date_discharged - LeftHemicolectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                  "90 day mortality" = mean(is.finite(date_death_ons) & date_death_ons - LeftHemicolectomy_date_admitted <= 90),
                                                  "90 day COVID-19" = mean(is.finite(LeftHemicolectomy_date) & LeftHemicolectomy_date - LeftHemicolectomy_date_admitted <= 90  & LeftHemicolectomy_date - LeftHemicolectomy_date_admitted >=0),
                                                 "90 day VTE" = mean(LeftHemicolectomypost.VTE, na.rm = T))],
      dt[is.finite(RightHemicolectomy_date_admitted)  & admit.wave.RightHemicolectomy == x,.("Procedures" = .N,
                                                       "Patients" = length(unique(patient_id)),
                                                       "Male" =  mean(sex=='M'),
                                                       "Age (IQR)" =  paste(quantile((RightHemicolectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                       "BMI (IQR)" =  paste(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                       "IMD (IQR)" =  paste(quantile(imd,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                       "1st Vaccination" = mean(is.finite(covid_vaccine_dates_1) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_1 >= 14),
                                                       "2nd Vaccination" = mean(is.finite(covid_vaccine_dates_2) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_2 >= 14),
                                                       "3rd Vaccination" = mean(is.finite(covid_vaccine_dates_3) & RightHemicolectomy_date_admitted  - covid_vaccine_dates_3 >= 14),
                                                       "Emergency" =  mean(substr(RightHemicolectomy_admission_method,1,1) == "2"),
                                                       "Length of stay (IQR)" =  paste(quantile((RightHemicolectomy_date_discharged - RightHemicolectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                       "90 day mortality" = mean(is.finite(date_death_ons) & date_death_ons - RightHemicolectomy_date_discharged <= 90),
                                                       "90 day COVID-19" = mean(is.finite(RightHemicolectomy_date) & RightHemicolectomy_date - RightHemicolectomy_date_admitted <= 90 & RightHemicolectomy_date - RightHemicolectomy_date_admitted >= 0),
                                                       "90 day VTE" = mean(RightHemicolectomypost.VTE, na.rm = T))],
      dt[is.finite(TotalColectomy_date_admitted)  & admit.wave.TotalColectomy == x,.("Procedures" = .N,
                                                   "Patients" =   length(unique(patient_id)),
                                                   "Male" =   mean(sex=='M'),
                                                   "Age (IQR)" =    paste(quantile((TotalColectomy_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                   "BMI (IQR)" =    paste(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                   "IMD (IQR)" =    paste(quantile(imd,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                   "1st Vaccination" = mean(is.finite(covid_vaccine_dates_1) & TotalColectomy_date_admitted  - covid_vaccine_dates_1 >= 14),
                                                   "2nd Vaccination" = mean(is.finite(covid_vaccine_dates_2) & TotalColectomy_date_admitted  - covid_vaccine_dates_2 >= 14),
                                                   "3rd Vaccination" = mean(is.finite(covid_vaccine_dates_3) & TotalColectomy_date_admitted  - covid_vaccine_dates_3 >= 14),
                                                   "Emergency" =    mean(substr(TotalColectomy_admission_method,1,1) == "2"),
                                                   "Length of stay (IQR)" =  paste(quantile((TotalColectomy_date_discharged - TotalColectomy_date_admitted),c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                   "90 day mortality" = mean(is.finite(date_death_ons) & date_death_ons - TotalColectomy_date_admitted <= 90),
                                                  "90 day COVID-19" = mean(is.finite(TotalColectomy_date) & TotalColectomy_date - TotalColectomy_date_admitted <= 90  & TotalColectomy_date - TotalColectomy_date_admitted >=0),
                                                  "90 day VTE" = mean(TotalColectomypost.VTE, na.rm = T))],
      dt[is.finite(RectalResection_date_admitted)  & admit.wave.RectalResection == x,.("Procedures" = .N,
                                                    "Patients" =   length(unique(patient_id)),
                                                    "Male" =   mean(sex=='M'),                                                      
                                                    "Age (IQR)" =   paste(quantile((RectalResection_date_admitted - dateofbirth)/365.25,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                    "BMI (IQR)" =   paste(quantile(bmi,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                    "IMD (IQR)" =   paste(quantile(imd,c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                    "1st Vaccination" = mean(is.finite(covid_vaccine_dates_1) & RectalResection_date_admitted  - covid_vaccine_dates_1 >= 14),
                                                    "2nd Vaccination" = mean(is.finite(covid_vaccine_dates_2) & RectalResection_date_admitted  - covid_vaccine_dates_2 >= 14),
                                                    "3rd Vaccination" = mean(is.finite(covid_vaccine_dates_3) & RectalResection_date_admitted  - covid_vaccine_dates_3 >= 14),
                                                    "Emergency" =   mean(substr(RectalResection_admission_method,1,1) == "2"),
                                                    "Length of stay (IQR)" =  paste(quantile((RectalResection_date_discharged - RectalResection_date_admitted),c(0.25,0.5,0.75),na.rm = T),collapse = ","),
                                                   "90 day mortality" = mean(is.finite(date_death_ons) & date_death_ons - RectalResection_date_admitted <= 90),
                                                   "90 day COVID-19" = mean(is.finite(RectalResection_date) & RectalResection_date - RectalResection_date_admitted <= 90  & RectalResection_date - RectalResection_date_admitted >=0),
                                                   "90 day VTE" = mean(RectalResectionpost.VTE, na.rm = T))])),
      cbind(dt[is.finite(LeftHemicolectomy_date_admitted) & admit.wave.LeftHemicolectomy == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))],
            dt[is.finite(RightHemicolectomy_date_admitted) & admit.wave.RightHemicolectomy == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))],
            dt[is.finite(TotalColectomy_date_admitted) & admit.wave.TotalColectomy == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))],
            dt[is.finite(RectalResection_date_admitted) & admit.wave.RectalResection == x,.N, keyby = region][, do.call(paste,c(.SD, sep = ": "))]),
cbind(dt[admit.wave.LeftHemicolectomy == x,.N, by = LeftHemicolectomy_primary_diagnosis][order(-N), do.call(paste,c(.SD, sep = ": "))][1:10],
      dt[admit.wave.RightHemicolectomy == x,.N, by = RightHemicolectomy_primary_diagnosis][order(-N), do.call(paste,c(.SD, sep = ": "))][1:10],
      dt[admit.wave.TotalColectomy == x,.N, by = TotalColectomy_primary_diagnosis][order(-N), do.call(paste,c(.SD, sep = ": "))][1:10],
      dt[admit.wave.RectalResection == x,.N, by = RectalResection_primary_diagnosis][order(-N), do.call(paste,c(.SD, sep = ": "))][1:10])), digits = 2), type = 'html', here::here("output",paste0("table1",x,".html"))))
  

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

dt[,(paste0(procedures,"_end_fu")) := lapply(.SD, function(x) data.table::fifelse(is.finite(date_death_ons) & x+90 > date_death_ons, date_death_ons,x+90)),
   .SDcols = paste0(procedures, '_date_discharged')]

dt[,no.op := rowSums(!is.na(.SD)),.SDcols = c(paste(procedures,"_date_admitted",sep = ""))]

dt[,max.date := lapply(.SD,max, na.rm = T),
   by = patient_id,
   .SDcols = c(paste0(procedures,"_end_fu"))]

dt[is.finite(gp.end) & max.date > gp.end, max.date := gp.end]
dt[!is.finite(max.date), max.date := data.table::as.IDate('2022-02-01')]


dt.fixed <- dt[,.SD, .SDcols = c(fixed,'max.date')]
dt.tv <- survival::tmerge(dt.fixed,dt.fixed,id = patient_id, end = event(max.date) )

dt.times <- dt[,.SD, .SDcols = c('patient_id',time.cols,"gp.end", proc.time.cols,paste(procedures,"_end_fu",sep = ""))]
dt.times[, gp.end := as.numeric(gp.end)]
dt.tv.values <- dt[,.SD, .SDcols = c('patient_id',
                                     paste(procedures,"_date_admitted",sep = ""),
                                     paste(procedures,"_emergency_readmit_date_admitted",sep = ""), 
                                     proc.tval.cols)]

for (i in c(proc.time.cols,paste0(procedures,"_end_fu"),"gp.end"))  {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.times,",
                           "id = patient_id,",
                           i," = event(",i,",",i,")))")))
}


for (i in time.cols) {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.times,",
                           "id = patient_id,",
                           i," = tdc(",i,",",i,",NA)))")))
}


for (proc in procedures) {
  for (val in proc.tval.stubs) {
    time.var = '_date_admitted'
    if (val == '_emergency_readmit_primary_diagnosis') time.var = '_emergency_readmit_date_admitted'
       eval(parse(text = paste0("assign(x = 'dt.tv', value = survival::tmerge(dt.tv,dt.tv.values,",
                             "id = patient_id,", 
                             paste0(proc,val)," = tdc(",paste0(proc,time.var),",",paste0(proc,val),",NA)))")))
  }
}

data.table::setDT(dt.tv)
for (i in c(proc.time.cols,paste0(procedures,"_end_fu")))  {
  eval(parse(text = paste0("assign(x = 'dt.tv', value = dt.tv[",i,"==0, ",i,":=NA])")))
}

dt.tv[, gp.end := max(gp.end, na.rm = T), by = patient_id]
dt.tv[gp.end == 0, gp.end := Inf]
data.table::setkey(dt.tv,patient_id, tstart, tstop)



admission.dates <- c('admit.date','discharge.date','end.fu')
dt.tv[,admit.date := do.call(pmax, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_admitted")]
dt.tv[!is.finite(admit.date), admit.date := NA]
dt.tv[,end.fu := do.call(pmin, c(.SD, na.rm = T)), .SDcols = c(paste0(procedures,"_end_fu"),"gp.end")] ## gp.end
dt.tv[!is.finite(end.fu), end.fu := NA]
dt.tv[,discharge.date := do.call(pmax, c(.SD, na.rm = T)), .SDcols = paste0(procedures,"_date_discharged")]
dt.tv[!is.finite(discharge.date), discharge.date := NA]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, (c('discharge.date','end.fu')) := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = c('discharge.date','end.fu')]
dt.tv[discharge.date > end.fu, discharge.date := NA]
dt.tv[,study.start := min(admit.date, na.rm = T), keyby = .(patient_id,end.fu)]
dt.tv[!is.finite(study.start), study.start := NA]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, (admission.dates) := lapply(.SD, data.table::nafill, type = "locf"), by = patient_id, .SDcols = admission.dates]

dt.tv[admit.date > discharge.date, (paste0(procedures,c('_admission_method','_primary_diagnosis',
                                                        '_days_in_critical_care',
                                                        '_case_category')))  := NA]

dt.tv[admit.date > discharge.date | is.na(admit.date), c('admit.date','discharge.date') := NA]

for (proc in procedures) {
  cols <-paste0(proc, proc.time.stubs)
  eval(parse(text = paste0("dt.tv[",proc,"_date_admitted>",proc,"_date_discharged,(cols) := NA]")))
}

dt.tv[!is.finite(admit.date), admit.date := NA]
dt.tv[!is.finite(end.fu), end.fu := NA]
dt.tv[!is.finite(discharge.date), discharge.date := end.fu]
dt.tv[!is.finite(study.start), study.start := NA]

## Assumption can't have two colonic resections on same day
dt.tv[admit.date == LeftHemicolectomy_date_admitted,op.type := 'LeftHemicolectomy']
dt.tv[admit.date == RightHemicolectomy_date_admitted,op.type := 'RightHemicolectomy']
dt.tv[admit.date == TotalColectomy_date_admitted,op.type := 'TotalColectomy']
dt.tv[admit.date == RectalResection_date_admitted,op.type := 'RectalResection']


data.table::setkey(dt.tv,patient_id,tstart,tstop)
id_change = dt.tv[, c(TRUE, patient_id[-1] != patient_id[-.N] | end.fu[-1] != end.fu[-.N])]
dt.tv[, op.type := lapply(.SD, function(x) x[cummax(((!is.na(x)) | id_change) * .I)]), .SDcols = "op.type"]

#dt.tv[!is.finite(op.type), op.type := NA]


data.table::setkey(dt.tv,patient_id,tstart,tstop)


dt.tv[,`:=`(admission_method = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_admission_method")]
dt.tv[,`:=`(primary_diagnosis = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_primary_diagnosis")]
dt.tv[,`:=`(days_in_critical_care = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_days_in_critical_care")]
dt.tv[,`:=`(case_category = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_case_category")]
dt.tv[,`:=`(emergency_readmit_primary_diagnosis = data.table::fcoalesce(.SD)), .SDcols = paste0(procedures,"_emergency_readmit_primary_diagnosis")]
dt.tv[,Current.Cancer := substr(primary_diagnosis,1,1) =='C']
dt.tv[is.na(Current.Cancer), Current.Cancer := F]
dt.tv[,(proc.tval.cols) := NULL]

dt.tv[,COVIDpositivedate  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(COVIDpositivedate), COVIDpositivedate := NA]
dt.tv[,emergency_readmitdate  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_emergency_readmit_date_admitted"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(emergency_readmitdate), emergency_readmitdate := NA]
dt.tv[,VTE_HES  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_VTE_HES_date_admitted"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(VTE_HES), VTE_HES := NA]
dt.tv[,VTE_GP  :=  min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_VTE_GP_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(VTE_GP), VTE_GP := NA]
dt.tv[,anticoagulation_prescriptions  := min(do.call(pmax, c(.SD, na.rm = T)), na.rm = T), .SDcols = paste0(procedures,"_anticoagulation_prescriptions_date"), by = .(patient_id, end.fu)]
dt.tv[!is.finite(anticoagulation_prescriptions), anticoagulation_prescriptions := NA]

dt.tv[,(proc.time.cols) := NULL]


###############################
# Pre operative risk factors----
##############################
dt.tv[,age := floor((tstart - as.numeric(as.Date(paste0(dob,'-15'))))/365.25)]
dt.tv[,age.cat := cut(age, breaks = c(18,50,70,80,90),ordered_result = T , right = T, include.lowest = T)]


### Calculate Charlson index at time of operation - tdc so date present from first recording
comorb.cols <- c(names(dt)[grep("^pre",names(dt))])

data.table::setkey(dt.tv,"patient_id","tstart","tstop")

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

dt.tv[,Charl12 := cut(Charlson, breaks = c(0,1,2,100), labels = c("None","Single","Multiple or Severe"), ordered = T)]

## Operation type


## Define cancer operations
#dt[,Cancer.Surgery := !is.na(surgery_cancer)] ## TODO  make this more specific to cancer related to operation type
#dt[is.na(Cancer.Surgery), Cancer.Surgery := F]

### Define elective or emergency operations
dt.tv[,Emergency := substr(as.character(admission_method),1,1)=="2"]
dt.tv[is.na(Emergency), Emergency := F]

## Define vaccination status - 14 days post date as effective
dt.tv[, vaccination.status := is.finite(covid_vaccine_dates_1) + is.finite(covid_vaccine_dates_2) + is.finite(covid_vaccine_dates_3)]


##############################
#Post operative outcomes----
##############################

### Length of stay----
dt.tv[,discharged := is.finite(data.table::shift(discharge.date, n = 1L, type = 'lead')) & data.table::shift(discharge.date, n = 1L, type = 'lead') == tstop]

### Post operative VTE----
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

dt.tv[, anticoagulation_prescriptions := lapply(.SD, data.table::nafill, type = "nocb"), by = patient_id, .SDcols = 'anticoagulation_prescriptions']

dt.tv[, post.VTE := ((is.finite(data.table::shift(VTE_GP, n = 1L, type = 'lead')) &  
                        data.table::shift(VTE_GP, n = 1L, type = 'lead') == tstop) | 
                       (is.finite(data.table::shift(VTE_HES, n = 1L, type = 'lead')) & 
                          data.table::shift(VTE_HES, n = 1L, type = 'lead') == tstop)) & 
        is.finite(data.table::shift(anticoagulation_prescriptions, n = 1L, type = 'lead')) & 
        is.finite(anticoagulation_prescriptions) & 
        tstop <= end.fu]  # events flagged at end of episode
dt.tv[post.VTE == T, post.VTE.date := tstop]
dt.tv[,post.VTE.date := min(post.VTE.date, na.rm = T), by = patient_id]


### Post operative Covid-19----
dt.tv[,COVIDpositive := is.finite(data.table::shift(COVIDpositivedate, n = 1L, type = 'lead')) & data.table::shift(COVIDpositivedate, n = 1L, type = 'lead') == tstop]
dt.tv[,postcovid := cumsum(COVIDpositive), by = .(patient_id, end.fu)]

### Readmissions----
dt.tv[,emergency_readmit  := is.finite(data.table::shift(emergency_readmitdate, n = 1L, type = 'lead')) & data.table::shift(emergency_readmitdate, n = 1L, type = 'lead') == tstop]

## Define types of emergency readmissions


### Mortality----
dt.tv[,died := tstop == date_death_ons]
## Cause of death


#######################
#Crude survival plots----
######################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

dt.cohort <- dt.tv[admit.date > end.fu | tstop > end.fu | tstart <admit.date | tstart < study.start]

dt.tv[, year := data.table::year(data.table::as.IDate(admit.date))]
dt.tv[, wave := factor(data.table::fifelse(admit.date <= as.numeric(data.table::as.IDate("2020-09-01")),"Wave_1",
                                    data.table::fifelse(admit.date <= as.numeric(data.table::as.IDate("2021-05-01")),"Wave_2",
                                                        data.table::fifelse(admit.date <= as.numeric(data.table::as.IDate("2021-12-31")),"Wave_3","Wave_4"))), ordered = F)]

dt.tv[, `:=`(start = tstart - study.start,
             end = tstop - study.start)]

###Counts
dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,.N,by = op.type]
dt.tv[died == 1,tail(.SD,1), by = .(patient_id,op.type)][,.N,by = op.type]
dt.tv[COVIDpositive == 1,tail(.SD,1), by = .(patient_id,op.type)][,.N,by = op.type]
dt.tv[post.VTE == 1,tail(.SD,1), by = .(patient_id,op.type)][,.N,by = op.type]
dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,mean(discharge.date - admit.date, na.rm = T),by = op.type]
dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,mean(Charlson, na.rm = T),by = op.type]

###

crude.surv <- survival::survfit(survival::Surv(start,end,died) ~ op.type, data = dt.tv[start>=0], id = patient_id)
plot_surv <- survminer::ggsurvplot(crude.surv, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_surv$plot,filename = "plot_surv.png", path=here::here("output"))

crude.covid <- survival::survfit(survival::Surv(start,end,COVIDpositive) ~ op.type, data = dt.tv[start>=0], id = patient_id)
plot_covid <- survminer::ggsurvplot(crude.covid, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_covid$plot,filename = "plot_covid.png", path=here::here("output"))

crude.los <- survival::survfit(survival::Surv(start,end,discharged) ~ op.type, data = dt.tv[is.finite(admit.date) & start>=0], id = patient_id)
plot_los <-survminer::ggsurvplot(crude.los, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_los$plot,filename = "plot_los.png", path=here::here("output"))


dt.tv[start >= 0,discharge.start := min(discharge.date, na.rm = T), by = .(patient_id, end.fu)]

crude.readmit <- survival::survfit(survival::Surv(tstart - discharge.start,tstop - discharge.start ,emergency_readmit) ~ op.type, data = dt.tv[tstart - discharge.start >=0], id = patient_id)
plot_readmit <-survminer::ggsurvplot(crude.readmit, data = dt.tv, risk.table = T)
ggplot2::ggsave(plot = plot_readmit$plot,filename = "plot_readmit.png", path=here::here("output"))

##############################
# Post operative COVID risk
#############################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

##TODO need to see model output to determine appropriate model fit.

post.op.covid.model <- survival::coxph(survival::Surv(start,end,COVIDpositive) ~ op.type + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer  + Emergency + Charl12, id = patient_id, data = dt.tv[start>=0 & year < 2022])
data.table::fwrite(broom::tidy(post.op.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_covid_model.csv"))



################################
# COVID impact on post operative Mortality
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.post.covid.covid.model <- survival::coxph(survival::Surv(start,end,died) ~ op.type + postcovid + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, id = patient_id, data = dt.tv[start>=0 & wave != 'Wave_4'])
data.table::fwrite(broom::tidy(post.op.post.covid.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_post_covid_model.csv"))


################################
# COVID impact on LOS
#################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.los.post.covid.model <- survival::coxph(survival::Surv(start,end,discharged) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, id = patient_id, data = dt.tv[start>=0 & !is.na(admit.date) & wave != 'Wave_4'])
data.table::fwrite(broom::tidy(post.op.los.post.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post.op.los.post.covid.model.csv"))

