library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T)
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
dt.tv[, sub.op := (is.finite(Colectomy) & Colectomy ==T) |
        (is.finite(Cholecystectomy) & Cholecystectomy == T) |
        (is.finite(HipReplacement)  & HipReplacement == T) | 
        (is.finite(KneeReplacement) & KneeReplacement == T)]

data.table::setkey(dt.tv,patient_id,tstart,tstop)
max.grp.col_(dt = 'dt.tv',
             max.var.name = 'sub.op',
             aggregate.cols = 'sub.op',
             id.vars = c("patient_id","end.fu"))

data.table::setkey(dt.tv,patient_id,tstart)

n.type.events <- sort(unique(dt.tv[(postop.readmit.cohort) & sub.op == T,event.readmit]))[-1]


post.op.readmit.model.sub <- 
  lapply(n.type.events, function(i) survival::coxph(survival::Surv(start.readmit,end.readmit,event.readmit==i) ~ Colectomy*wave + Cholecystectomy*wave + 
                                            KneeReplacement*wave + postcovid*wave + sex +  age.cat + bmi.cat + imd5 + vaccination.status.factor + Current.Cancer + 
                                              Emergency + LOS.bin + Charl12 + recentCOVID + previousCOVID + region, id = patient_id,
                                          data = dt.tv[(postop.readmit.cohort) & sub.op == T], model = T))

#data.table::fwrite(broom::tidy(post.op.readmit.model.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopreadmitmodelsub.csv"))
names(post.op.readmit.model.sub) <- c('COVID-19','Non COVID-19 emergency readmission','Mortality')[n.type.events]
modelsummary::modelsummary(post.op.readmit.model.sub,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, conf_level = .95, exponentiate = TRUE, output = here::here("output","postopreadmitmodelsub.html"))


new.data.postop.covid <- data.table::data.table('start.readmit' = rep(0,8*length(procedures)),
                                                'end.readmit' = rep(30,8*length(procedures)),
                                                'event.readmit' = rep(F,8*length(procedures)),
                                                'Colectomy' = c(rep(T,8),rep(F,24)),
                                                'Cholecystectomy'=c(rep(F,8),rep(T,8),rep(F,16)),
                                                'HipReplacement'=c(rep(F,16),rep(T,8),rep(F,8)),
                                                'KneeReplacement'=c(rep(F,24),rep(T,8)),
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

cuminc.adjusted.readmit.sub <- 
  matrix(cuminc.cox(n.type.events = n.type.events,dt = 'dt.tv[(postop.readmit.cohort) & sub.op == T]', model = 'post.op.readmit.model.sub', newdata = 'new.data.postop.covid', day = 90), byrow = T, ncol = 4)

colnames(cuminc.adjusted.readmit.sub) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.readmit.sub) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.readmit.model.sub,cuminc.adjusted.readmit.sub, file = here::here("output","postopreadmit_sub.RData"))
data.table::fwrite(cuminc.adjusted.readmit.sub, file = here::here("output","postopreadmit_sub.csv"))


readmit.waves.sub.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.readmit.sub, keep.rownames = T),
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

ggplot2::ggsave(plot = readmit.waves.sub.plot, here::here('output','readmit_waves_sub_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )