library(foreach)
library(data.table)
library(modelsummary)

ncores <- parallel::detectCores(logical = T) 
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

###########################################################

dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))
procedures <- c('Colectomy','Cholecystectomy',
                'HipReplacement','KneeReplacement')


data.table::setkey(dt.tv,patient_id,tstart,tstop)


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

post.op.died.model.sub <- 
  list(survival::coxph(survival::Surv(start,end,died) ~ Colectomy*wave + Cholecystectomy*wave + KneeReplacement*wave + 
                         postcovid*wave +
                         Emergency*wave + 
                          age.cat + sex  + imd5 + vaccination.status.factor + Current.Cancer +
                         Charl12 + recentCOVID + previousCOVID, id = patient_id,
                       data = dt.tv[start >=0 & sub.op == T & any.op == T], model = T))
#data.table::fwrite(broom::tidy(post.op.died.model.sub[[1]], exponentiate= T, conf.int = T), file = here::here("output","postopdiedmodelsub.csv"))


names(post.op.died.model.sub) <- c('Mortality')
modelsummary::modelsummary(post.op.died.model.sub,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopdiedmodelsub.html"))
modelsummary::modelsummary(post.op.died.model.sub,estimate  = "{estimate} [{conf.low}, {conf.high}], (p = {p.value})", statistic = NULL, gof_omit = "Num.Obs.|n|nevent", conf_level = .95, exponentiate = TRUE, output = here::here("output","postopdiedmodelsub.txt"))

  new.data.postop.covid <- data.table::data.table('start' = rep(0,16*length(procedures)),
                                                  'end' = rep(30,16*length(procedures)),
                                                  'event' = rep(F,16*length(procedures)),
                                                  'Colectomy' = c(rep(T,16),rep(F,48)),
                                                  'Cholecystectomy'=c(rep(F,16),rep(T,16),rep(F,32)),
                                                'HipReplacement'=c(rep(F,32),rep(T,16),rep(F,16)),
                                                'KneeReplacement'=c(rep(F,48),rep(T,16)),
                                                'postcovid' = rep(rep(c(rep(F,4),rep(T,4)),
                                                                  times = length(procedures)),each = 2),
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
cuminc.adjusted.mortality.sub.long <- cuminc.cox(n.type.events = n.type.events,
           dt = 'dt.tv[start >=0  & sub.op == T & any.op == T]',
           model = 'post.op.died.model.sub',
           newdata = 'new.data.postop.covid', day = 90)
cuminc.adjusted.mortality.sub <-   matrix(cuminc.adjusted.mortality.sub.long, byrow = T, ncol = 4)


colnames(cuminc.adjusted.mortality.sub) <- paste0('Wave_',1:4)
rownames(cuminc.adjusted.mortality.sub) <- paste0(rep(rep(c(' No COVID post ',' COVID post '),each = 2),times = 4),paste0(rep(c('Elective ','Emergency '), times = 8),rep(procedures, each = 4)))

save(post.op.died.model.sub,
cuminc.adjusted.mortality.sub,
cuminc.adjusted.mortality.sub.long,
new.data.postop.covid, 
file = here::here("output","postopmortality_sub.RData"))
data.table::fwrite(cuminc.adjusted.mortality.sub, file = here::here("output","postopmortality_sub.csv"))


cuminc.long.labelled.sub <-cbind(new.data.postop.covid,cuminc.adjusted.mortality.sub.long, c(rep('Colectomy',16), rep('Cholecystectomy',16), rep('HipReplacement',16), rep('KneeReplacement',16))) 
cuminc.long.labelled.sub[,wave := readr::parse_number(wave)]
cuminc.long.labelled.sub[,Emergency := factor(Emergency, levels = c(F,T), labels = c('Elective','Emergency'))]
names(cuminc.long.labelled.sub)[names(cuminc.long.labelled.sub)=='V3'] <- 'Operation'


mortality.waves.sub.plot <- ggplot2::ggplot(cuminc.long.labelled.sub[order(Emergency, Operation, postcovid,wave),], #data.table::melt(data.table::data.table(cuminc.adjusted.mortality.sub, keep.rownames = T),
                                                    #     id.vars = 'rn',
                                                     #    variable.name = 'Wave',
                                                      #   value.name = '90 Day Cumulative Mortality Risk (%)')[, `:=`(COVID = grepl('^COVID*',rn),
                                                       #                                                          Operation = gsub('^No COVID|^COVID', '',rn))],
                                            ggplot2::aes(x = wave, 
                                                         y = cuminc.adjusted.mortality.sub.long, 
                                                         linetype = postcovid,
                                                         colour = Operation
                                            )) + 
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~Emergency, scales = 'free') +
  ggplot2::ylab('90 Day Cumulative Mortality Risk (%)') +
  ggplot2::xlab('COVID-19 wave') 

ggplot2::ggsave(plot = mortality.waves.sub.plot, here::here('output','mortality_waves_sub_plot.png'),dpi = 'retina', width = 7, height = 5, units = 'in', device = 'png' )
