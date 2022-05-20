load(file = here::here("output","cohort_long.RData"))
procedures <- unique(dt.tv$op.type)

#######################
#Crude survival plots----
######################


###Counts for sense check----
data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,.N,by = op.type]
data.table::dcast(dt.tv[,tail(.SD,1), by = .(patient_id, op.type)], op.type ~ died, length)
data.table::dcast(dt.tv[,tail(.SD,1), by = .(patient_id, op.type)], op.type ~ postcovid, length)
data.table::dcast(dt.tv[,tail(.SD,1), by = .(patient_id, op.type)], op.type ~ postVTEany, length)
dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,mean(discharge.date - admit.date, na.rm = T),by = op.type]
dt.tv[,tail(.SD,1), by = .(patient_id,op.type)][,mean(Charlson, na.rm = T),by = op.type]


### 90 day table km estimates ----

covariates <- c('age.cat','sex','op.type','postcovid','Charl12','bmi.cat','imd5','region','vaccination.status.factor','Emergency','Current.Cancer','wave')

# Overall survival----
crude.surv.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[start>=0 ,as.factor(get(covariates[i]))]))),
                                                                           levels(dt.tv[start>=0 ,as.factor(get(covariates[i]))]),
                                                                           round(as.data.frame(summary(survival::survfit(survival::Surv(start,end,died) ~ get(covariates[i]),
                                                                                                                   data = dt.tv[start>=0], 
                                                                                                                   id = patient_id), 
                                                                                                 times = 90,
                                                                                                 extend = T)[c('n.risk','n.event','surv','lower','upper')]), digits = 3))))
crude.surv.cov <- crude.surv.cov[,`:=`(cuminc = 1 - surv,
                                         lower95 = 1 - upper,
                                         upper95 = 1 - lower)]
data.table::fwrite(crude.surv.cov,file = here::here("output", "surv_90day_counts.csv"))

# Covid risk----
crude.covid.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[start>=0 & tstop <= covid.end,as.factor(get(covariates[i]))]))),
                                                                           levels(dt.tv[start>=0 & tstop <= covid.end,as.factor(get(covariates[i]))]),
                                                                           round(as.data.frame(summary(survival::survfit(survival::Surv(start,end,COVIDpositive) ~ get(covariates[i]),
                                                                                                                   data = dt.tv[start>=0 & tstop <= covid.end], 
                                                                                                                   id = patient_id), 
                                                                                                 times = 90,
                                                                                                 extend = T)[c('n.risk','n.event','surv','lower','upper')]), digits = 3))))
crude.covid.cov <- crude.covid.cov[,`:=`(cuminc = 1 - surv,
                                         lower95 = 1 - upper,
                                         upper95 = 1 - lower)]
data.table::fwrite(crude.covid.cov,file = here::here("output", "covid_90day_counts.csv"))

# Readmit risk----
crude.readmit.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end,as.factor(get(covariates[i]))]))),
                                                                            levels(dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end,as.factor(get(covariates[i]))]),
                                                                            round(as.data.frame(summary(survival::survfit(survival::Surv(tstart - discharge.start,tstop - discharge.start ,emergency_readmit) ~ get(covariates[i]),
                                                                                                                    data = dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end], 
                                                                                                                    id = patient_id), 
                                                                                                  times = 90,
                                                                                                  extend = T)[c('n.risk','n.event','surv','lower','upper')]), digits = 3))))
crude.readmit.cov <- crude.readmit.cov[,`:=`(cuminc = 1 - surv,
                                         lower95 = 1 - upper,
                                         upper95 = 1 - lower)]
data.table::fwrite(crude.readmit.cov,file = here::here("output", "readmit_90day_counts.csv"))

# VTE risk----
crude.VTE.cov <- data.table::rbindlist(lapply(1:length(covariates), function(i) cbind(rep(covariates[i],length(levels(dt.tv[start >=0 & tstop <= VTE.end,as.factor(get(covariates[i]))]))),
                                                                              levels(dt.tv[start >=0 & tstop <= VTE.end,as.factor(get(covariates[i]))]),
                                                                              round(as.data.frame(summary(survival::survfit(survival::Surv(start, end ,post.VTE) ~ get(covariates[i]),
                                                                                                                      data = dt.tv[start >=0 & tstop <= VTE.end], 
                                                                                                                      id = patient_id), 
                                                                                                    times = 90,
                                                                                                    extend = T)[c('n.risk','n.event','surv','lower','upper')]), digits = 3))))
crude.VTE.cov <- crude.VTE.cov[,`:=`(cuminc = 1 - surv,
                                             lower95 = 1 - upper,
                                             upper95 = 1 - lower)]
data.table::fwrite(crude.VTE.cov,file = here::here("output", "VTE_90day_counts.csv"))


### Plots----
# By procedure----
# Overall survival----

crude.surv <- survival::survfit(survival::Surv(start,end,died) ~ op.type, data = dt.tv[start>=0], id = patient_id)
plot_surv <- survminer::ggsurvplot(crude.surv, data = dt.tv, risk.table = T,break.time.by = 10,  xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_surv),
                filename = "plot_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')

crude.covid <- survival::survfit(survival::Surv(start,end,COVIDpositive) ~ op.type, data = dt.tv[start>=0 & tstop <= covid.end], id = patient_id)
plot_covid <- survminer::ggsurvplot(crude.covid, data = dt.tv, risk.table = T,break.time.by = 10, xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_covid),filename = "plot_covid_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')

crude.los <- survival::survfit(survival::Surv(start,end,discharged) ~ op.type, data = dt.tv[is.finite(admit.date) & start>=0 & tstop <= los.end], id = patient_id)
plot_los <-survminer::ggsurvplot(crude.los, data = dt.tv, risk.table = T,break.time.by = 10, xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_los), filename = "plot_los_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')


crude.readmit <- survival::survfit(survival::Surv(tstart - discharge.start,tstop - discharge.start ,emergency_readmit) ~ op.type, data = dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end], id = patient_id)
plot_readmit <-survminer::ggsurvplot(crude.readmit, data = dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end], risk.table = T,break.time.by = 10, xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_readmit),filename = "plot_readmit_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')
# By wave

crude.surv.wave <- survival::survfit(survival::Surv(start,end,died) ~ wave, data = dt.tv[start>=0], id = patient_id)
plot_surv.wave <- survminer::ggsurvplot(crude.surv.wave, data = dt.tv, risk.table = T,break.time.by = 10, xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_surv.wave),
                filename = "plot_surv_wave_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')

crude.covid <- survival::survfit(survival::Surv(start,end,COVIDpositive) ~ wave, data = dt.tv[start>=0 & tstop <= covid.end], id = patient_id)
plot_covid.wave <- survminer::ggsurvplot(crude.covid, data = dt.tv, risk.table = T,break.time.by = 10,xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_covid.wave),filename = "plot_covid_wave_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')

crude.los.wave <- survival::survfit(survival::Surv(start,end,discharged) ~ wave, data = dt.tv[is.finite(admit.date) & start>=0 & tstop <= los.end], id = patient_id)
plot_los.wave <-survminer::ggsurvplot(crude.los.wave, data = dt.tv, risk.table = T,break.time.by = 10,xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_los.wave), filename = "plot_los_wave_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')


dt.tv[start >= 0,discharge.start := min(discharge.date, na.rm = T), by = .(patient_id, end.fu)]

crude.readmit.wave <- survival::survfit(survival::Surv(tstart - discharge.start,tstop - discharge.start ,emergency_readmit) ~ wave, data = dt.tv[tstart - discharge.start >=0 & tstop <=readmit.end], id = patient_id)
plot_readmit.wave <-survminer::ggsurvplot(crude.readmit.wave, data = dt.tv, risk.table = T,break.time.by = 10,xlim = c(0,90))
ggplot2::ggsave(plot = survminer:::.build_ggsurvplot(plot_readmit.wave),filename = "plot_readmit_wave_km.png", path=here::here("output"),
                dpi = 'retina', width = 7, height = 7, units = 'in')
