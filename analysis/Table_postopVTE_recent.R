library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

## start = 0 = operation / admit date , postcovid.VTE.cohort starts on days 1
data.table::setkey(dt.tv,patient_id,tstart,tstop)

# Not enough deaths to treat separately from emergency readmissions
dt.tv[event.VTE == 3, event.VTE := 2]
n.type.events <- sort(unique(dt.tv[(postcovid.VTE.cohort) ,event.VTE]))[-1]

post.op.VTE.model.recentCOVID <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.VTE==i) ~ Abdominal*wave + Obstetrics*wave + CardioThoracicVascular*wave + postcovid*wave   + age.cat +
                                            sex  + imd5 +  vaccination.status.factor + Current.Cancer + Emergency + LOS.bin + Charl12 + recentCOVID*wave + previousCOVID , id = patient_id,
                                          data = dt.tv[(postcovid.VTE.cohort)], model = T))

#data.table::fwrite(broom::tidy(post.op.VTE.model.recentCOVID[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopVTEmodelrecentCOVID.csv"))

names(post.op.VTE.model.recentCOVID) <- c('Post discharge VTE','Non COVID-19 emergency readmission or mortality')[n.type.events]
modelsummary::modelsummary(post.op.VTE.model.recentCOVID,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodelrecentCOVID.html"))
modelsummary::modelsummary(post.op.VTE.model.recentCOVID,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodelrecentCOVID.txt"))


new.data.postop.recent.covid <- data.table::data.table('start.readmit' = rep(0,8*length(procedures)),
                                                'end.readmit' = rep(90,8*length(procedures)),
                                                'event.VTE' = rep(F,8*length(procedures)),
                                                'Abdominal' = c(rep(T,8),rep(F,24)),
                                            #    'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
                                                'Obstetrics'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'Orthopaedic'=c(rep(F,16),rep(T,8),rep(F,8)),
                                             #   'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
                                                'CardioThoracicVascular'=c(rep(F,24),rep(T,8)),
                                                'postcovid' = rep(F,8*length(procedures)),
                                                'age.cat' = rep('(50,70]',8*length(procedures)),
                                                'sex' = rep('F',8*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                'vaccination.status.factor' = rep('3',8*length(procedures)),
                                                'region' = rep("East Midlands",8*length(procedures)),
                                                'Current.Cancer' = rep(T,8*length(procedures)),
                                                'Emergency' =   rep(F,8*length(procedures)),
                                                'LOS.bin' =   rep(F,8*length(procedures)),
                                                'Charl12' =  rep('Single',8*length(procedures)),
                                                'recentCOVID' =rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.recent.VTE <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort)]', 
                    model = 'post.op.VTE.model.recentCOVID', newdata = 'new.data.postop.recent.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.recent.VTE) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.recent.VTE) <- paste0(c('No recent COVID','Recent COVID'),rep(procedures, each = 2))

save(post.op.VTE.model.recentCOVID,cuminc.adjusted.recent.VTE, file = here::here("output","postopVTErecent.RData"))
data.table::fwrite(cuminc.adjusted.recent.VTE, file = here::here("output","postopVTE_recentCOVID.csv"))

VTE.waves.recent.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.recent.VTE, keep.rownames = T),
                                 id.vars = 'rn',
                                 variable.name = 'Wave',
                                 value.name = '90 Day Cumulative VTE Incidence (%)')[, `:=`(`Recent COVID` = grepl('^Recent COVID*',rn),
                                                                                            Operation = gsub('^No recent COVID|^Recent COVID', '',rn))],
                ggplot2::aes(x = Wave, 
                             y = `90 Day Cumulative VTE Incidence (%)`, 
                             group = rn,
                             colour = Operation,
                             linetype = `Recent COVID`)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = VTE.waves.recent.plot, here::here('output','VTE_waves_recentCOVID_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

