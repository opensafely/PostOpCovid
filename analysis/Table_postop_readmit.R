library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Orthopaedic',
'Obstetrics','CardioThoracicVascular')
dt.tv[(postop.readmit.cohort),.N ,keyby = .(wave, Abdominal, Obstetrics, Orthopaedic, CardioThoracicVascular, event.readmit)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)


for(v in procedures) {
  max.grp.col_(dt = 'dt.tv',
             max.var.name = v,
             aggregate.cols = v,
             id.vars = c("patient_id","end.fu"))
}
#dt.tv <- dt.tv[Abdominal == T | Obstetrics == T | Orthopaedic == T,]

data.table::setkey(dt.tv, patient_id, tstart)
n.type.events <- sort(unique(dt.tv[(postop.readmit.cohort) ,event.readmit]))[-1]




post.op.readmit.model <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.readmit==i) ~ Abdominal*wave + 
                                                   #  Cardiac*wave +
                                                     Obstetrics*wave  +
                                                    #   Thoracic*wave  + 
                                                     CardioThoracicVascular*wave  +
                                                      postcovid*wave +  
                                                       sex + age.cat + 
                                                      bmi.cat + imd5 +
                                                      vaccination.status.factor  + 
                                                       Current.Cancer + Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID + region, 
                                                    id = patient_id,
                                          data = dt.tv[(postop.readmit.cohort)], model = T))

names(post.op.readmit.model) <- c('Emergency Readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.readmit.model,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopreadmitmodel.html"))


#data.table::fwrite(broom::tidy(post.op.readmit.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopreadmitmodel1.csv"))
#data.table::fwrite(broom::glance(post.op.readmit.model[[1]]), file = here::here("output","postopreadmitmodelsummary1.csv"))
#data.table::fwrite(broom::tidy(post.op.readmit.model[[2]], exponentiate= T, conf.int = T), file = here::here("output","postopreadmitmodel2.csv"))
#data.table::fwrite(broom::glance(post.op.readmit.model[[2]]), file = here::here("output","postopreadmitmodelsummary2.csv"))


new.data.postop.covid <- data.table::data.table('start.readmit' = rep(0,8*length(procedures)),
                                                'end.readmit' = rep(30,8*length(procedures)),
                                                'event.readmit' = rep(F,8*length(procedures)),
                                                'Abdominal' = c(rep(T,8),rep(F,8*length(procedures) - 8)),
                                            #    'Cardiac'=c(rep(F,8*length(procedures) -5*8),rep(T,8),rep(F,8*length(procedures) -2*8)),
                                                'Obstetrics'=c(rep(F,8*length(procedures) -3*8),rep(T,8),rep(F,8*length(procedures) -2*8)),
                                                'Orthopaedic'=c(rep(F,8*length(procedures) -2*8),rep(T,8),rep(F,8*length(procedures) -3*8)),
                                             #  'Thoracic'=c(rep(F,8*length(procedures) - 2*8),rep(T,8),rep(F,8*length(procedures) -5*8)),
                                               'CardioThoracicVascular'=c(rep(F,8*length(procedures) - 8),rep(T,8)),
                                                'postcovid' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
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
                                                'recentCOVID' = rep(F,8*length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.readmit <- 
  matrix(cuminc.cox(n.type.events = n.type.events,
  dt = 'dt.tv[(postop.readmit.cohort)]',
   model = 'post.op.readmit.model', 
   newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.readmit) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.readmit) <-  paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.readmit.model,cuminc.adjusted.readmit, file = here::here("output","postopreadmit.RData"))
data.table::fwrite(cuminc.adjusted.readmit, file = here::here("output","postopreadmit.csv"))


readmit.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.readmit, keep.rownames = T),
                                 id.vars = 'rn',
                                 variable.name = 'Wave',
                                 value.name = '90 Day Cumulative Readmission Incidence (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                                Operation = gsub('^No COVID|^COVID', '',rn))],
                ggplot2::aes(x = Wave, 
                             y = `90 Day Cumulative Readmission Incidence (%)`, 
                             group = rn,
                             colour = Operation,
                             linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = readmit.waves.plot, here::here('output','readmit_waves_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

