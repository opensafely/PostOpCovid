library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Abdominal','Obstetrics','Orthopaedic','CardioThoracicVascular')

data.table::setkey(dt.tv,patient_id,tstart,tstop)

post.op.died.model <- 
  list(survival::coxph(survival::Surv(start,end,died) ~ Abdomina  + Obstetrics + CardioThoracicVascular + 
                         postcovid*wave + Emergency*wave + age.cat + sex + imd5  + vaccination.status.factor + Current.Cancer +
                        Charl12 + recentCOVID + previousCOVID, id = patient_id,
                       data = dt.tv[start >=0 & any.op == T], model = T))
#data.table::fwrite(broom::tidy(post.op.died.model[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopdiedmodel.csv"))

names(post.op.died.model) <- c('Mortality')
modelsummary::modelsummary(post.op.died.model,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopdiedmodel.html"))


new.data.postop.covid <- data.table::data.table('start' = rep(0,16*length(procedures)),
                                                'end' = rep(30,16*length(procedures)),
                                                'event' = rep(F,16*length(procedures)),
                                                'Abdominal' = c(rep(T,16),rep(F,48)),
                                          #      'Cardiac'=c(rep(F,8),rep(T,8),rep(F,32)),
                                                'Obstetrics'=c(rep(F,16),rep(T,16),rep(F,32)),
                                                'Orthopaedic'=c(rep(F,32),rep(T,16),rep(F,16)),
                                           #     'Thoracic'=c(rep(F,32),rep(T,8),rep(F,8)),
                                                'CardioThoracicVascular'=c(rep(F,48),rep(T,16)),
                                                'postcovid' = rep(rep(c(rep(F,4),rep(T,4)), times = length(procedures)),each = 2),
                                                'age.cat' = rep('(50,70]',16*length(procedures)),
                                                'sex' = rep('F',16*length(procedures)),
                                                'bmi.cat' = rep(levels(dt.tv$bmi.cat)[2],16*length(procedures)),
                                                'imd5' = rep(levels(dt.tv$imd5)[3], 16*length(procedures)),
                                                'wave' = rep(paste0('Wave_',1:4),times = 4*length(procedures)),
                                                'vaccination.status.factor' = rep('3',16*length(procedures)),
                                                'region' = rep("East Midlands",16*length(procedures)),
                                                'Current.Cancer' = rep(T,16*length(procedures)),
                                                'Emergency' =   rep(rep(c(F,T), each = 4),2*length(procedures)),
                                                'LOS.bin' =   rep(F,16*length(procedures)),
                                                'Charl12' =  rep('Single',16*length(procedures)),
                                                'recentCOVID' = rep(F,16*length(procedures)),
                                                'previousCOVID' = rep(F,16*length(procedures)),
                                                'patient_id' = 1:(16*length(procedures)))

n.type.events <- 1

cuminc.adjusted.mortality.long <- cuminc.cox(n.type.events = n.type.events,
           dt = 'dt.tv[start >=0 & any.op == T]',
           model = 'post.op.died.model',
           newdata = 'new.data.postop.covid', 
           day = 90)

cuminc.adjusted.mortality <-   matrix(cuminc.adjusted.mortality.long, byrow = T, ncol = 4)


colnames(cuminc.adjusted.mortality) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.mortality) <- paste0(rep(rep(c(' No COVID post ',' COVID post '),each = 2),times = 4),paste0(rep(c('Elective ','Emergency '), times = 8),rep(procedures, each = 4)))

save(post.op.died.model,cuminc.adjusted.mortality, file = here::here("output","postopmortality.RData"))
data.table::fwrite(cuminc.adjusted.mortality, file = here::here("output","postopmortality.csv"))

cuminc.long.labelled <-cbind(new.data.postop.covid,cuminc.adjusted.mortality.long, c(rep('Abdominal',16), rep('Obstetrics',16), rep('Orthopaedic',16), rep('CardioThoracicVascular',16))) 
cuminc.long.labelled[,wave := readr::parse_number(wave)]
cuminc.long.labelled[,Emergency := factor(Emergency, levels = c(F,T), labels = c('Elective','Emergency'))]
names(cuminc.long.labelled)[names(cuminc.long.labelled)=='V3'] <- 'Operation'


mortality.waves.plot <- ggplot2::ggplot(cuminc.long.labelled[order(Emergency, Operation, postcovid,wave)], #data.table::melt(data.table::data.table(cuminc.adjusted.mortality, keep.rownames = T),
                                 #id.vars = 'rn',
                                 #variable.name = 'Wave',
                                 #value.name = '90 Day Cumulative Mortality Risk (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                  #                                                       Operation = gsub('^No COVID|^COVID', '',rn))],
                ggplot2::aes(x = wave, 
                             y = cuminc.adjusted.mortality.long, 
                             linetype = postcovid,
                             colour = Operation
                              )) + 
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~Emergency, scales = 'free') +
  ggplot2::ylab('90 Day Cumulative Mortality Risk (%)') +
  ggplot2::xlab('COVID-19 wave') 

ggplot2::ggsave(plot = mortality.waves.plot, here::here('output','mortality_waves_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )

