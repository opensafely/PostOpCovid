library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular', 'Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement','Major.op','Intermediate.or.Minor.op')

## Count variables for demographic tables
dt.tv[,postVTE90.perepisode := event.VTE == 1 & (postcovid.VTE.cohort)]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'postVTE90.perepisode',aggregate.cols = 'postVTE90.perepisode',id.vars = c("patient_id", "end.fu"))
dt.tv[!is.finite(postVTE90.perepisode),postVTE90.perepisode := 0]

dt.tv[,postCOVID30.perepisode := event == 1 &  (postop.covid.cohort) & end <= 30]
data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',max.var.name = 'postCOVID30.perepisode',aggregate.cols = 'postCOVID30.perepisode',id.vars = c("patient_id", "end.fu"))
dt.tv[!is.finite(postCOVID30.perepisode),postCOVID30.perepisode := 0]
dt.tv[date_death_ons - discharge.date <=90, death_underlying_cause_ons := NA]

n.ops <- rnd(dt.tv[start ==0  & is.finite(admit.date) & any.op == T, lapply(.SD,function(x) sum(x == T)), .SDcols = c(procedures)])

n.pats <- rnd(length(unique(dt.tv[start ==0 & is.finite(admit.date) & any.op == T,patient_id])))

start.date <- dt.tv[start ==0 & is.finite(admit.date) & any.op == T, min(as.Date(as.integer(admit.date), origin = as.Date('1970-01-01')))]

last.date <-  dt.tv[start ==0  & is.finite(admit.date) & any.op == T, max(as.Date(as.integer(admit.date), origin = as.Date('1970-01-01')))]

demo.tab <- 
  data.table::transpose(
    cbind(
    data.table::data.table("procedures" = procedures),
                              foreach::foreach(i = 1:length(procedures), .combine = 'rbind', .inorder = T) %do% dt.tv[start ==0 & 
                              final.date >= tstart  & 
                              final.date >= tstop &
                               end <= 90 & 
                               any.op == T &
                                get(paste0(procedures[i])) == T,head(.SD,1),keyby = .(patient_id,end.fu)][,
                                                                                                                      .("Procedures" = rnd(.N),
                                                                                                                        "Patients" = rnd(length(unique(patient_id))),
                                                                                                                        "Female" = n.perc(sex=='F',dig = 3),
                                                                                                                        "Age (IQR)" = median.iqr(age,dig = 0),
                                                                                                                        "BMI < 19" =  n.perc(bmi.cat ==  '(1,18]', dig  = 3),
                                                                                                                        "BMI 19-24" =  n.perc(bmi.cat ==  '(18,24]', dig  = 3),
                                                                                                                        "BMI 25-29" =  n.perc(bmi.cat ==  '(24,29]', dig  = 3),
                                                                                                                        "BMI > 30" =  n.perc(bmi.cat ==  '(29,100]', dig  = 3),
                                                                                                                        "BMI missing" = n.perc(bmi.cat ==  'Missing', dig  = 3),
                                                                                                                        "IMD quintile 1" = n.perc(as.numeric(imd5) ==  1, dig  = 3),
                                                                                                                        "IMD quintile 2" = n.perc(as.numeric(imd5) ==  2, dig  = 3),
                                                                                                                        "IMD quintile 3" = n.perc(as.numeric(imd5) ==  3, dig  = 3),
                                                                                                                        "IMD quintile 4" = n.perc(as.numeric(imd5) ==  4, dig  = 3),
                                                                                                                        "IMD quintile 5" = n.perc(as.numeric(imd5) ==  5, dig  = 3),
                                                                                                                        "IMD quintile Missing" = n.perc(imd5 == 'Missing', dig  = 3),
                                                                                                                        "Region: East" = n.perc(region ==  'East', dig  = 3),
                                                                                                                        "Region: East Midlands" = n.perc(region ==  'East Midlands', dig  = 3),
                                                                                                                        "Region: London" = n.perc(region ==  'London', dig  = 3),
                                                                                                                        "Region: North East" = n.perc(region ==  'North East', dig  = 3),
                                                                                                                        "Region: North West" = n.perc(region ==  'North West', dig  = 3),
                                                                                                                        "Region: South East" = n.perc(region == 'South East', dig  = 3),
                                                                                                                        "Region: South West" = n.perc(region == 'South West', dig  = 3),
                                                                                                                        "Region: West Midlands" = n.perc(region == 'West Midlands', dig  = 3),
                                                                                                                        "Region: Yorkshire and The Humber" = n.perc(region == 'Yorkshire and The Humber', dig  = 3),
                                                                                                                        "Wave 1" = n.perc(wave == 'Wave_1', dig = 3),
                                                                                                                        "Wave 2" = n.perc(wave == 'Wave_2', dig = 3),
                                                                                                                        "Wave 3" = n.perc(wave == 'Wave_3', dig = 3),
                                                                                                                        "Wave 4" = n.perc(wave == 'Wave_4', dig = 3),
                                                                                                                        "1st Vaccination" = n.perc(vaccination.status.factor==1,dig = 3),
                                                                                                                        "2nd Vaccination" = n.perc(vaccination.status.factor==2,dig = 3),
                                                                                                                        "3rd Vaccination" = n.perc(vaccination.status.factor==3,dig = 3),
                                                                                                                        "Current Cancer"  = n.perc(substr(primary_diagnosis,1,1) =='C',dig = 3),
                                                                                                                        "Emergency" = n.perc(substr(admission_method,1,1) == "2",dig = 3),
                                                                                                                        "Charlson index" = median.iqr(Charlson, dig = 1),
                                                                                                                        "Length of stay (IQR)" =  median.iqr(discharge.date - admit.date,dig = 0),
                                                                                                                        "90 day emergency readmission (%)" = n.perc(date.90.day(x = emergency_readmitdate, ref.dat = discharge.date + 1), dig = 3),
                                                                                                                        "90 day mortality (%)" =n.perc(date.90.day(x = date_death_ons, ref.dat = admit.date),dig = 4),
                                                                                                                        "30 day post operative COVID-19 (%)" = n.perc(postCOVID30.perepisode, dig = 3),
                                                                                                                        "Critical care on index admission (%)" = n.perc(is.finite(days_in_critical_care) & days_in_critical_care > 0, dig = 3),
                                                                                                                        "Recent 7 to 42 days COVID-19 (%)" = n.perc(recentCOVID, dig = 3),
                                                                                                                        "Previous > 42 COVID-19 (%)" = n.perc(previousCOVID, dig = 3),
                                                                                                                        "90 day VTE post discharge (%)" = n.perc(postVTE90.perepisode, dig = 4))]
  #                             data.table::transpose(
  #                               cbind(paste0("Primary admission diagnosis category - ",1:5),
  #                                     dt.tv[(postop.covid.cohort) & start ==0,
  #                                           lapply(.SD,function(x) rnd(sum(x))),
  #                                           keyby = c('Prim.admit.cat'), .SDcols = c(procedures)][,
  #                                                                                                    lapply(.SD,
  #                                                                                                           function(x) paste0(Prim.admit.cat,
  #                                                                                                                              ": ",
  #                                                                                                                              x,
  #                                                                                                                              ' (',
  #                                                                                                                              round(100*x/sum(x,na.rm = T),
  #                                                                                                                                    digits = 1),
  #                                                                                                                              '%)')[order(-x)]),
  #                                                                                                    .SDcols = 2:(length(procedures) + 1)][1:5,]),
  #                               make.names = 'V1'),
  #                             data.table::transpose(
  #                               cbind(paste0("Emergency readmission diagnosis category excluding COVID - ",1:5),
  #                                     dt.tv[(postop.readmit.cohort)  & start.readmit > 0 &
  #                                           event.readmit == 1 & !is.na(Readmit.cat) & COVIDreadmission != T,
  #                                           lapply(.SD,function(x) rnd(sum(x))),
  #                                           keyby = c('Readmit.cat'), .SDcols = c(procedures)][,
  #                                                                                                    lapply(.SD,
  #                                                                                                           function(x) paste0(Readmit.cat,
  #                                                                                                                              ": ",
  #                                                                                                                              x,
  #                                                                                                                              ' (',
  #                                                                                                                              round(100*x/sum(x,na.rm = T),
  #                                                                                                                                    digits = 1),
  #                                                                                                                              '%)')[order(-x)]),
  #                                                                                                    .SDcols = 2:(length(procedures) + 1)][1:5,]),
  #                               make.names = 'V1'),
  #                             data.table::transpose(
  #                               cbind(paste0("Underlying cause of death category - ",1:5),
  #                                     dt.tv[end.fu == tstop & !is.na(COD.cat) ,
  #                                           lapply(.SD,function(x) rnd(sum(x))),
  #                                           keyby = c('COD.cat'), .SDcols = c(procedures)][,
  #                                                                                                                      lapply(.SD,
  #                                                                                                                             function(x) paste0(COD.cat,
  #                                                                                                                                                ": ",
  #                                                                                                                                                x,
  #                                                                                                                                                ' (',
  #                                                                                                                                                round(100*x/sum(x,na.rm = T),
  #                                                                                                                                                      digits = 1),
  #                                                                                                                                                '%)')[order(-x)]),
  #                                                                                                                      .SDcols = 2:(length(procedures) + 1)][1:5,]),
  #                               make.names = 'V1')
   ),
  make.names = 'procedures',
  keep.names = 'Characteristics')

demo.tab[, ]


data.table::fwrite(demo.tab, file = here::here("output","table_demo.csv"))


