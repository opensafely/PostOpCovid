load(file = here::here("output","cohort_long.RData"))
procedures <- unique(dt.tv$op.type)
################################
# Post operative COVID risk ----
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")
covariates <- c('age.cat','sex','op.type','postcovid','Charl12','bmi.cat','imd5','region','vaccination.status.factor','Emergency','Current.Cancer','wave','recentCOVID','previousCOVID')

post.op.covid.overall.model <- 
  survival::coxph(survival::Surv(start,end,COVIDpositive) ~ wave + age + sex + bmi.cat + imd5 + region + vaccination.status.factor + Current.Cancer + Emergency + Charlson + recentCOVID + previousCOVID, id = patient_id,
                  data = dt.tv[start>=0 & tstop <= covid.end  ], model = T)
data.table::fwrite(broom::tidy(post.op.covid.overall.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_covid_overall_model.csv"))


post.op.covid.model <- 
  survival::coxph(survival::Surv(start,end,COVIDpositive) ~ op.type + wave + age + sex + bmi.cat + imd5 + region + vaccination.status.factor + Current.Cancer + Emergency + Charlson + recentCOVID + previousCOVID, id = patient_id,
                  data = dt.tv[start>=0 & tstop <= covid.end  ], model = T)
data.table::fwrite(broom::tidy(post.op.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_covid_model.csv"))

covid.risk.30day <- predict(object = post.op.covid.model, 
                                   newdata = data.table::data.table('start' = rep(0,8),
                                                                    'end' = rep(30,8*length(procedures)),
                                                                    'COVIDpositive' = rep(F,8*length(procedures)),
                                                                    'recentCOVID' = rep(0,8*length(procedures)),
                                                                    'previousCOVID' = rep(0,8*length(procedures)),
                                                                    'op.type' = rep(procedures,each = 8),
                                                                    'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                                    'age' = rep(60,8*length(procedures)),
                                                                    'sex' = rep('M',8*length(procedures)),
                                                                    'region' = rep("East Midlands",8*length(procedures)),
                                                                    'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                                    'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                                    'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                                    'Current.Cancer' = rep(T,8*length(procedures)),
                                                                    'Emergency' =  rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                                    'Charlson' =  rep(1,8*length(procedures)),
                                                                    'patient_id' = 1:8*length(procedures)), type = 'expected',se.fit = T)
covid.risk.ci.30day <- matrix(paste0(round((1- exp(-covid.risk.30day$fit))*100,3),
                                   ' (', round((1 - exp(-(covid.risk.30day$fit - 1.96*covid.risk.30day$se.fit)))*100,3),',',
                                   round((1 - exp(-(covid.risk.30day$fit + 1.96*covid.risk.30day$se.fit)))*100,3),')'),nrow = 4)

rownames(covid.risk.ci.30day) <- paste0('Wave_',1:4)
colnames(covid.risk.ci.30day) <- paste0(c('Elective_','Emergency_'),rep(procedures, each = 2))

covid.risk.ci.30day

data.table::fwrite(covid.risk.ci.30day,file = here::here("output", "post_op_covid_cuminc_model.csv"))

#flexmodelcovid <- flexsurv::flexsurvreg(survival::Surv(start,end,COVIDpositive) ~ op.type + wave + age + sex + bmi + vaccination.status.factor + Current.Cancer + Emergency + Charlson,
#                     data = dt.tv[start>=0 & tstop <= covid.end  ], model = T, dist = 'gengamma')


#plot(flexmodelcovid)


print(xtable::xtable(finalfit::finalfit.coxph(dt.tv[start>=0 & tstop <= covid.end  ],
  'survival::Surv(start,end,COVIDpositive)',
   covariates
)),type = 'html',here::here("output","post_op_covid_ff_model.html"))



################################
# Post operative VTE risk----
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.VTE.model <- 
  survival::coxph(survival::Surv(start,end,post.VTE) ~ op.type + wave + postcovid + age + sex + bmi.cat + imd5  + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                  data = dt.tv[start>=0 & tstop <= VTE.end  ])
data.table::fwrite(broom::tidy(post.op.VTE.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_VTE_model.csv"))

VTE.risk.30day <- predict(object = post.op.VTE.model, 
                                   newdata = data.table::data.table('start' = rep(0,8),
                                                                    'end' = rep(30,8*length(procedures)),
                                                                    'post.VTE' = rep(F,8*length(procedures)),
                                                                    'postcovid' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                                    'recentCOVID' = rep(0,8*length(procedures)),
                                                                    'previousCOVID' = rep(0,8*length(procedures)),
                                                                    'op.type' = rep(procedures,each = 8),
                                                                    'vaccination.status' = rep(3,8*length(procedures)),
                                                                    'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                                    'age' = rep(60,8*length(procedures)),
                                                                    'sex' = rep('M',8*length(procedures)),
                                                                    'region' = rep("East Midlands",8*length(procedures)),
                                                                    'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                                    'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                                    'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                                    'Current.Cancer' = rep(T,8*length(procedures)),
                                                                    'Emergency' =   rep(F,8*length(procedures)), 
                                                                    'Charl12' =  rep('Single',8*length(procedures)),
                                                                    'patient_id' = 1:8*length(procedures)), type = 'expected', se.fit = T)
VTE.risk.ci.30day <- matrix(paste0(round((1- exp(-VTE.risk.30day$fit))*100,3),
                                   ' (', round((1- exp(-(VTE.risk.30day$fit - 1.96*VTE.risk.30day$se.fit)))*100,3),',',
                                               round((1- exp(-(VTE.risk.30day$fit + 1.96*VTE.risk.30day$se.fit)))*100,3),')'),nrow = 4)
  
rownames(VTE.risk.ci.30day) <- paste0('Wave_',1:4)
colnames(VTE.risk.ci.30day) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

VTE.risk.ci.30day

data.table::fwrite(VTE.risk.ci.30day,file = here::here("output", "post_op_VTE_cuminc_model.csv"))


print(xtable::xtable(finalfit::finalfit.coxph(dt.tv[start>=0 & tstop <= VTE.end],
                                                        'survival::Surv(start,end,post.VTE)',
                                                        covariates
)),type = 'html',here::here("output","post_op_VTE_ff_model.html"))

################################
# COVID impact on post operative Mortality ----
##################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.post.covid.surv.model <- 
  survival::coxph(survival::Surv(start,end,died) ~ op.type + postcovid + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                  data = dt.tv[start>=0 ])
data.table::fwrite(broom::tidy(post.op.post.covid.surv.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_post_covid_surv_model.csv"))



post.op.post.covid.surv.waves.model <- 
  survival::coxph(survival::Surv(start,end,died) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                  data = dt.tv[start>=0 ])
data.table::fwrite(broom::tidy(post.op.post.covid.surv.waves.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_post_covid_surv_waves_model.csv"))


print(xtable::xtable(finalfit::finalfit.coxph(dt.tv[start>=0 ],
                                                        'survival::Surv(start,end,died)',
                                                        covariates
)), type = 'html',here::here("output","post_op_mort_ff_model.html"))

################################
# COVID impact on LOS ----
#################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.los.post.covid.model <- survival::coxph(survival::Surv(start,end,discharged) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, 
                                                id = patient_id, data = dt.tv[start>=0 & tstop <= los.end & !is.na(admit.date) ])
data.table::fwrite(broom::tidy(post.op.los.post.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_los_post_covid_model.csv"))

print(xtable::xtable(finalfit::finalfit.coxph(dt.tv[start>=0 & tstop <= los.end & !is.na(admit.date) ],
                                                        'survival::Surv(start,end,discharged)',
                                                        covariates
)), type = 'html', file = here::here("output","post_op_los_ff_model.html"))


################################
# COVID impact on Emergency Readmission ----
#################################
data.table::setkey(dt.tv,"patient_id","tstart","tstop")

post.op.los.post.covid.model <- survival::coxph(survival::Surv(start,end,emergency_readmit) ~ op.type + postcovid*wave + age + sex + bmi + factor(vaccination.status, ordered = F) + Current.Cancer + Emergency + Charl12 + recentCOVID + previousCOVID, 
                                                id = patient_id, data = dt.tv[start>=0 & tstop <= readmit.end])
data.table::fwrite(broom::tidy(post.op.los.post.covid.model, exponentiate= T, conf.int = T), file = here::here("output","post_op_readmit_post_covid_model.csv"))

print(xtable::xtable(finalfit::finalfit.coxph(dt.tv[start>=0 & tstop <=  readmit.end ],
                                                        'survival::Surv(start,end,emergency_readmit)',
                                                        covariates
)), type = 'html', file = here::here("output","post_op_readmit_ff_model.html"))