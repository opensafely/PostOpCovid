load(file = here::here("output","cohort_long.RData"))
procedures <- unique(dt.tv$op.type)
################################
# Post operative COVID risk
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.covid.model <- 
  survival::coxph(survival::Surv(start,end,COVIDpositive) ~ op.type + wave + age + sex + bmi + vaccination.status.factor + Current.Cancer + Emergency + Charlson, id = patient_id,
                  data = dt.tv[start>=0 & tstop <= covid.end  ], model = T)
data.table::fwrite(broom::tidy(post.op.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_covid_model.csv"))

dt.tv[,.(op.type,wave,age,sex,bmi,vaccination.status.factor,Current.Cancer,Emergency,Charl12)]

covid.risk.30day <- matrix(predict(object = post.op.covid.model, 
                                   newdata = data.table::data.table('start' = rep(0,8),
                                                                    'end' = rep(30,8*length(procedures)),
                                                                    'COVIDpositive' = rep(F,8*length(procedures)),
                                                                    'op.type' = rep(procedures,each = 8),
                                                                    'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                                    'age' = rep(60,8*length(procedures)),
                                                                    'sex' = rep('M',8*length(procedures)),
                                                                    'bmi' = rep(25,8*length(procedures)),
                                                                    'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                                    'Current.Cancer' = rep(T,8*length(procedures)),
                                                                    'Emergency' =  rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                                    'Charlson' =  rep(1,8*length(procedures)),
                                                                    'patient_id' = 1:8*length(procedures)), type = 'expected'), nrow = 4)

rownames(covid.risk.30day) <- paste0('Wave_',1:4)
colnames(covid.risk.30day) <- paste0(c('Elective_','Emergency_'),rep(procedures, each = 2))

data.table::fwrite(round(covid.risk.30day*100,1),file = here::here("output", "post_op_covid_cuminc_model.csv"))

#flexmodelcovid <- flexsurv::flexsurvreg(survival::Surv(start,end,COVIDpositive) ~ op.type + wave + age + sex + bmi + vaccination.status.factor + Current.Cancer + Emergency + Charlson,
#                     data = dt.tv[start>=0 & tstop <= covid.end  ], model = T, dist = 'gengamma')


#plot(flexmodelcovid)



################################
# Post operative VTE risk
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.VTE.model <- 
  survival::coxph(survival::Surv(start,end,post.VTE) ~ op.type + wave + postcovid + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, id = patient_id,
                  data = dt.tv[start>=0 & tstop <= VTE.end  ])
data.table::fwrite(broom::tidy(post.op.VTE.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_VTE_model.csv"))


################################
# COVID impact on post operative Mortality
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.post.covid.covid.model <- 
  survival::coxph(survival::Surv(start,end,died) ~ op.type + postcovid + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, id = patient_id,
                  data = dt.tv[start>=0 ])
data.table::fwrite(broom::tidy(post.op.post.covid.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_post_covid_model.csv"))



post.op.post.covid.covid.waves.model <- 
  survival::coxph(survival::Surv(start,end,died) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, id = patient_id,
                  data = dt.tv[start>=0 ])
data.table::fwrite(broom::tidy(post.op.post.covid.covid.waves.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_post_covid_waves_model.csv"))



################################
# COVID impact on LOS
#################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.los.post.covid.model <- survival::coxph(survival::Surv(start,end,discharged) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12, 
                                                id = patient_id, data = dt.tv[start>=0 & tstop <= los.end & !is.na(admit.date) ])
data.table::fwrite(broom::tidy(post.op.los.post.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_los_post_covid_model.csv"))

