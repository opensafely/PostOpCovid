library(foreach)
library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

data.table::setkey(dt.tv,patient_id,tstart,tstop)

post.op.died.model <- 
  list(survival::coxph(survival::Surv(start,end,died) ~ Abdominal + Cardiac + Obstetrics + Thoracic + Vascular + 
                         postcovid*wave + age.cat + sex + bmi.cat + imd5  + vaccination.status.factor + region + Current.Cancer +
                         Emergency + Charl12 + recentCOVID + previousCOVID, id = patient_id,
                       data = dt.tv[start >=0 & any.op == T], model = T))
data.table::fwrite(broom::tidy(post.op.died.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopdiedmodel.csv"))


new.data.postop.covid <- data.table::data.table('start' = rep(0,8*length(procedures)),
                                                'end' = rep(30,8*length(procedures)),
                                                'event' = rep(F,8*length(procedures)),
                                                'Abdominal' = c(rep(T,8),rep(F,40)),
                                                'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
                                                'Obstetrics'=c(rep(F,16),rep(T,8),rep(F,24)),
                                                'Orthopaedic'=c(rep(F,24),rep(T,8),rep(F,16)),
                                                'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
                                                'Vascular'=c(rep(F,40),rep(T,8)),
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

n.type.events <- 1
cuminc.adjusted.mortality <-   matrix(cuminc.cox(n.type.events = n.type.events,
                                                 dt = 'dt.tv[start >=0 & any.op == T]',
                                                 model = 'post.op.died.model',
                                                 newdata = 'new.data.postop.covid', 
                                                 day = 90), byrow = T, ncol = 4)


colnames(cuminc.adjusted.mortality) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.mortality) <- paste0(c('No COVID','COVID'),rep(procedures, each = 2))

save(post.op.died.model,cuminc.adjusted.mortality, file = here::here("output","postopmortality.RData"))
data.table::fwrite(cuminc.adjusted.mortality, file = here::here("output","postopmortality.csv"))

mortality.waves.plot <- ggplot2::ggplot(data.table::melt(data.table::data.table(cuminc.adjusted.mortality, keep.rownames = T),
                                 id.vars = 'rn',
                                 variable.name = 'Wave',
                                 value.name = '90 Day Cumulative Mortality Risk (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                                                         Operation = gsub('^No COVID|^COVID', '',rn))],
                ggplot2::aes(x = Wave, 
                             y = `90 Day Cumulative Mortality Risk (%)`, 
                             group = rn,
                             colour = Operation,
                             linetype = COVID)) +
  ggplot2::geom_line()

ggplot2::ggsave(plot = mortality.waves.plot, here::here('output','mortality_waves_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

