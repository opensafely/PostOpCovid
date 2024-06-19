library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')


data.table::setkey(dt.tv,patient_id,tstart,tstop)

dt.tv[, sub.op :=( (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T))& Emergency== F]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))

# Not enough deaths to treat separately from emergency readmissions
dt.tv[event.VTE == 3, event.VTE := 2]
n.type.events <- sort(unique(dt.tv[(postcovid.VTE.cohort)  & sub.op == T,event.VTE]))[-1]

post.op.VTE.model <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start,end,event.VTE==i) ~ Colectomy*wave + Cholecystectomy*wave +  HipReplacement*wave + KneeReplacement*wave + postcovid*wave  + age.cat +
                                                      sex +  imd5 + vaccination.status.factor  + Current.Cancer + Charl12 + recentCOVID*wave, id = patient_id,
                                                    data = dt.tv[(postcovid.VTE.cohort) & sub.op == T], model = T))

#data.table::fwrite(broom::tidy(post.op.VTE.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopVTEmodel.csv"))

names(post.op.VTE.model) <- c('Post discharge VTE','Non COVID-19 emergency readmission or mortality')[n.type.events]
modelsummary::modelsummary(post.op.VTE.model,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", 
                           statistic = NULL, gof_omit = "Num.Obs.|n|nevent",
                           conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodel_procedure_specific.html"))
modelsummary::modelsummary(post.op.VTE.model,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", 
                           statistic = NULL, gof_omit = "Num.Obs.|n|nevent",
                           conf_level = .95, exponentiate = TRUE, output = here::here("output","postopVTEmodel_procedure_specific.txt"))
## with vaccination postcovid
new.data.postop.covid <- data.table::data.table('start.readmit' = rep(0,8*length(procedures)),
                                                'end.readmit' = rep(90,8*length(procedures)),
                                                'event.VTE' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
                                                'recentCOVID' =  rep(F,8*length(procedures)),
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
                                                'postcovid' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.VTE <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & sub.op == T]', model = 'post.op.VTE.model', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.VTE) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.VTE) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.VTE.model,cuminc.adjusted.VTE, file = here::here("output","postopVTEmodel_procedure_specific.RData"))
data.table::fwrite(cuminc.adjusted.VTE, file = here::here("output","postopVTEmodel_procedure_specific.csv"))

VTE.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.VTE, keep.rownames = T),
                                                   id.vars = 'rn',
                                                   variable.name = 'Wave',
                                                   value.name = '90 Day Cumulative VTE Incidence (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                                              Operation = gsub('^No COVID|^COVID', '',rn))],
                                  ggplot2::aes(x = Wave, 
                                               y = `90 Day Cumulative VTE Incidence (%)`, 
                                               group = rn,
                                               colour = Operation,
                                               linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = VTE.waves.plot, here::here('output','VTE_waves_procedure_specific.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

## without vaccination post covid
new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(90,8*length(procedures)),
                                                'event.VTE' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
                                                'recentCOVID' =  rep(F,8*length(procedures)),
                                                'age.cat' = rep('(50,70]',8*length(procedures)),
                                                'sex' = rep('F',8*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                'vaccination.status.factor' = rep('0',8*length(procedures)),
                                                'region' = rep("East Midlands",8*length(procedures)),
                                                'Current.Cancer' = rep(T,8*length(procedures)),
                                                'Emergency' =   rep(F,8*length(procedures)),
                                                'LOS.bin' =   rep(F,8*length(procedures)),
                                                'Charl12' =  rep('Single',8*length(procedures)),
                                                'postcovid' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.VTE <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & sub.op == T]', model = 'post.op.VTE.model', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.VTE) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.VTE) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.VTE.model,cuminc.adjusted.VTE, file = here::here("output","postopVTE_procedure_specific_unvac.RData"))
data.table::fwrite(cuminc.adjusted.VTE, file = here::here("output","postopVTE_procedure_specific_unvac.csv"))

VTE.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.VTE, keep.rownames = T),
                                                   id.vars = 'rn',
                                                   variable.name = 'Wave',
                                                   value.name = '90 Day Cumulative VTE Incidence (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                                              Operation = gsub('^No COVID|^COVID', '',rn))],
                                  ggplot2::aes(x = Wave, 
                                               y = `90 Day Cumulative VTE Incidence (%)`, 
                                               group = rn,
                                               colour = Operation,
                                               linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = VTE.waves.plot, here::here('output','VTE_waves_procedure_specific_unvac.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

## with vaccination recent covid
new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(90,8*length(procedures)),
                                                'event.VTE' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
                                                'postcovid' =  rep(F,8*length(procedures)),
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
                                                'recentCOVID' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.VTE <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & sub.op == T]', model = 'post.op.VTE.model', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.VTE) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.VTE) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.VTE.model,cuminc.adjusted.VTE, file = here::here("output","recentcovid_procedure_specific.RData"))
data.table::fwrite(cuminc.adjusted.VTE, file = here::here("output","recentcovid_procedure_specific.csv"))

VTE.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.VTE, keep.rownames = T),
                                                   id.vars = 'rn',
                                                   variable.name = 'Wave',
                                                   value.name = '90 Day Cumulative VTE Incidence (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                                              Operation = gsub('^No COVID|^COVID', '',rn))],
                                  ggplot2::aes(x = Wave, 
                                               y = `90 Day Cumulative VTE Incidence (%)`, 
                                               group = rn,
                                               colour = Operation,
                                               linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = VTE.waves.plot, here::here('output','recentcovid_VTE_procedure_specific.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

## without vaccination recent covid
new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(90,8*length(procedures)),
                                                'event.VTE' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
                                                'postcovid' =  rep(F,8*length(procedures)),
                                                'age.cat' = rep('(50,70]',8*length(procedures)),
                                                'sex' = rep('F',8*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],8*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 8*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 2*length(procedures)),
                                                'vaccination.status.factor' = rep('0',8*length(procedures)),
                                                'region' = rep("East Midlands",8*length(procedures)),
                                                'Current.Cancer' = rep(T,8*length(procedures)),
                                                'Emergency' =   rep(F,8*length(procedures)),
                                                'LOS.bin' =   rep(F,8*length(procedures)),
                                                'Charl12' =  rep('Single',8*length(procedures)),
                                                'recentCOVID' = rep(c(rep(F,4),rep(T,4)), times = length(procedures)),
                                                'previousCOVID' = rep(F,8*length(procedures)),
                                                'patient_id' = 1:(8*length(procedures)))

cuminc.adjusted.VTE <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postcovid.VTE.cohort) & sub.op == T]', model = 'post.op.VTE.model', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.VTE) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.VTE) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.VTE.model,cuminc.adjusted.VTE, file = here::here("output","recentcovid_postopVTE_procedure_specific_unvac.RData"))
data.table::fwrite(cuminc.adjusted.VTE, file = here::here("output","recentcovid_postopVTE_procedure_specific_unvac.csv"))

VTE.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.VTE, keep.rownames = T),
                                                   id.vars = 'rn',
                                                   variable.name = 'Wave',
                                                   value.name = '90 Day Cumulative VTE Incidence (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                                              Operation = gsub('^No COVID|^COVID', '',rn))],
                                  ggplot2::aes(x = Wave, 
                                               y = `90 Day Cumulative VTE Incidence (%)`, 
                                               group = rn,
                                               colour = Operation,
                                               linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = VTE.waves.plot, here::here('output','recentcovid_VTE_waves_procedure_specific_unvac.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

