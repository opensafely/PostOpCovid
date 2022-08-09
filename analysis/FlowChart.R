library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

procedures <- c('Abdominal','Cardiac','Obstetrics','Orthopaedic','Thoracic', 'Vascular')
source(here::here("analysis","Utils.R"))
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

data.table::setkey(dt.tv,patient_id,tstart,tstop)

##########################

n.ops <- rnd(dt.tv[is.finite(end.fu) & start ==0 & is.finite(admit.date) & any.op == T  & admit.date <= end.fu & final.date >= tstart,tail(.SD,1), keyby = .(patient_id, end.fu)][,lapply(.SD,function(x) sum(x == T)), .SDcols = c(procedures)])
n.ops.VTE <- rnd(dt.tv[(postcovid.VTE.cohort) & start.readmit >0   & end.readmit <=90, tail(.SD,1), keyby = .(patient_id, end.fu)][,lapply(.SD,function(x) sum(x == T)), .SDcols = c(procedures)])
n.ops.COVID <- rnd(dt.tv[(postop.covid.cohort) & start ==0  & final.date >= tstart , tail(.SD,1), keyby = .(patient_id, end.fu)][,lapply(.SD,function(x) sum(x == T)), .SDcols = c(procedures)])

n.pats <- rnd(length(unique(dt.tv[,patient_id])))
n.pats.study <- rnd(length(unique(dt.tv[is.finite(end.fu) &  start ==0 & is.finite(admit.date) & any.op == T & admit.date <= end.fu & final.date >= tstart,patient_id])))
n.pats.late <- rnd(length(unique(dt.tv[is.finite(end.fu) &  (admit.date > gp.end | admit.date >= end.fu | !is.finite(gp.end)) & (!is.finite(final.date) |   final.date <= tstart) & any.op == T,patient_id])))

n.covid.90 <- rnd(dt.tv[(postop.covid.cohort) ,max(event == 1 & start >=0 & end <=90 & any.op.COVID == T, na.rm = T) , keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.covid.90.censored <- rnd(dt.tv[(postop.covid.cohort),max(event == 1 & start >=0  & end <=90 & tstop <= final.date  & any.op.COVID == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.covid.30.censored <- rnd(dt.tv[(postop.covid.cohort),max(event == 1 & start >=0  & end <=30 & tstop <= final.date & any.op.COVID == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.covid.7.censored <- rnd(dt.tv[(postop.covid.cohort),max(event == 1 & start >=0  &  end <=7 & tstop <= final.date & any.op.COVID == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])

n.VTE.90 <- rnd(dt.tv[(postcovid.VTE.cohort),max(event.VTE == 1 & start.readmit >0   & end.readmit <=90 , na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.VTE.90.censored <- rnd(dt.tv[ (postcovid.VTE.cohort),max(event.VTE == 1 & start.readmit >0  & end.readmit <=90  & tstop <= final.date.VTE & any.op.VTE == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.VTE.30.censored <- rnd(dt.tv[(postcovid.VTE.cohort),max(event.VTE == 1 & start.readmit >0  & end.readmit <=30 & tstop <= final.date.VTE & any.op.VTE == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.VTE.7.censored <- rnd(dt.tv[(postcovid.VTE.cohort),max(event.VTE == 1 & start.readmit >0  & end.readmit <=7 & tstop <= final.date.VTE & any.op.VTE == T, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])

n.surv.90 <- rnd(dt.tv[start >= 0 & any.op == T & tstop <= end.fu,max(died == 1 & end <=90, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.surv.90.censored <- rnd(dt.tv[start >= 0 & any.op == T & tstop <= end.fu,max(died == 1 & end <=90, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.surv.30.censored <- rnd(dt.tv[start >= 0 & any.op == T & tstop <= end.fu,max(died == 1 & end <=30, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])
n.surv.7.censored <- rnd(dt.tv[start >= 0 & any.op == T  & tstop <= end.fu,max(died == 1 & end <=7, na.rm = T), keyby = .(patient_id, end.fu)][,tail(.SD,1), keyby = .(patient_id, end.fu)][,sum(V1==1)])


start.date <- dt.tv[ start ==0 & is.finite(admit.date) & any.op == T, min(as.Date(as.integer(admit.date), origin = as.Date('1970-01-01')))]

last.date <-  dt.tv[ start ==0  & is.finite(admit.date) & any.op == T, max(as.Date(as.integer(admit.date), origin = as.Date('1970-01-01')))]

p <- ggplot2::ggplot(data = data.frame(x = 1:100, y = 1:100), ggplot2::aes(x , y)) + 
  ggplot2::scale_x_continuous(minor_breaks = seq(10, 100, 10)) +
  ggplot2::scale_y_continuous(minor_breaks = seq(10, 100, 10)) +   
  ggplot2::theme_void()
  
p <- p + ggplot2::geom_rect(xmin = 32, xmax=67, ymin=94, ymax=100, color='black',
            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=97,label= paste0(n.pats,' Patients with coded procedures \n between ',start.date,' & ', last.date), size=2.5) 
  
p <- p + ggplot2::geom_rect(xmin = 32, xmax=68, ymin=81.5, ymax=90, color='black',
            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=86,label= paste0(n.pats.study,' Patients with procedures within \n registration period for a \n primary care practice contributing \n to SystmOne and no COVID 7 days prior'), size=2.5) 
  
p <- p + ggplot2::geom_rect(xmin = 70, xmax=101, ymin=86, ymax=98, color='black',
            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 85, y=92,label= paste0(n.pats - n.pats.study,' Patients excluded: \n Procedures when not registered \n to a primary care practice \n contributing to SystmOne'), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 28, xmax=71, ymin=62, ymax=79, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=71,label= paste0(sum(n.ops)," Procedures in categories \n(not mutually exclusive):\n",paste(paste0(names(n.ops),':'),n.ops, collapse = ' procedures \n '),' procedures'), size=2.5)

p <- p +
  ggplot2::geom_segment(
    x=50, xend=50, y=94, yend=90, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  ggplot2::geom_segment(
    x=50, xend=69.7, y=92, yend=92, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  ggplot2::geom_segment(
    x=50, xend=50, y=81.5, yend=79, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  ggplot2::geom_segment(
    #middle arrow 
    x=50, xend=50, y=62, yend=55, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Left arrow
  ggplot2::geom_segment(
    x=17, xend=17, y=58.5, yend=55, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Right arrow
  ggplot2::geom_segment(
    x=83, xend=83, y=58.5, yend=55, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Horizontal line
  ggplot2::geom_segment(
    x=17, xend=83, y=58.5, yend=58.5, 
    size=0.15, linejoin = "mitre", lineend = "butt")

##########################################
# COVID counts

p <- p + ggplot2::geom_rect(xmin = 0, xmax=33, ymin=49, ymax=55, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 17, y=52.2,label= paste0(sum(n.ops.COVID)," Procedures \n with no COVID-19 positivity \n in prior 7 days"), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 0, xmax=33, ymin=41, ymax=47, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 17, y=44.2,label= paste0("COVID-19 positivity within 90 days: \n ",n.covid.90," events\n",round(100 * n.covid.90/sum(n.ops.COVID), digits = 2),"%" ), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 0, xmax=33, ymin=28, ymax=39, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 17, y=33.2,label= paste0("COVID-19 positivity within 90 days \n censored for \n non COVID-19 readmissions: \n ",n.covid.90.censored," events\n",round(100 * n.covid.90.censored/sum(n.ops.COVID), digits = 2),"%" ), size=2.5)

p <- p + ggplot2::geom_rect(xmin =0, xmax=33, ymin=13, ymax=26, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 17, y=19.2,label= paste0("COVID-19 positivity within 30 days \n censored for \n non COVID-19 readmissions: \n ",n.covid.30.censored," events\n",round(100 * n.covid.30.censored/sum(n.ops.COVID), digits = 2),"%" ), size=2.5)

#######################################
# VTE counts

p <- p + ggplot2::geom_rect(xmin = 34, xmax=66, ymin=49, ymax=55, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=52.2,label= paste0(sum(n.ops.VTE)," Procedures \n for post discharge VTE \n analysis"), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 34, xmax=66, ymin=41, ymax=47, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=44.2,label= paste0("VTE 90 days of discharge: \n ",n.VTE.90," events\n",round(100 * n.VTE.90/sum(n.ops.VTE), digits = 2),"%" ), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 34, xmax=66, ymin=28, ymax=39, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=33.2,label= paste0("VTE 90 days of discharge\n censored for \n readmissions: \n ",n.VTE.90.censored," events\n",round(100 * n.VTE.90.censored/sum(n.ops.VTE), digits = 2),"%" ), size=2.5)

p <- p + ggplot2::geom_rect(xmin =34, xmax=66, ymin=13, ymax=26, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 50, y=19.2,label= paste0("VTE 30 days of discharge\n censored for \n readmissions:  \n ",n.VTE.30.censored," events\n",round(100 * n.VTE.30.censored/sum(n.ops.VTE), digits = 2),"%" ), size=2.5)

#########################################
# VTE counts

p <- p + ggplot2::geom_rect(xmin = 67, xmax=100, ymin=49, ymax=55, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 83, y=52.2,label= paste0(sum(n.ops)," Procedures \n for post COVID mortality\n analysis"), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 67, xmax=100, ymin=41, ymax=47, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 83, y=44.2,label= paste0("Deaths within 90 days: \n ",n.surv.90," events\n",round(100 * n.surv.90/sum(n.ops), digits = 2),"%" ), size=2.5)

p <- p + ggplot2::geom_rect(xmin = 67, xmax=100, ymin=28, ymax=39, color='black',
                            fill='white', size=0.25) +
  ggplot2::annotate('text', x= 83, y=33.2,label= paste0("Deaths within 30 days\n ",n.surv.30.censored," events\n",round(100 * n.surv.30.censored/sum(n.ops), digits = 2),"%" ), size=2.5)


########################################
p <- p +
  ggplot2::geom_segment(
    # Centre top arrow
    x=50, xend=50, y=49, yend=47, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Left top arrow
  ggplot2::geom_segment(
    x=17, xend=17, y=49, yend=47, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Right top arrow
  ggplot2::geom_segment(
    x=83, xend=83, y=49, yend=47, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  ggplot2::geom_segment(
    # Centre middle arrow
    x=50, xend=50, y=41, yend=39, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  #Left middle arrow
  ggplot2::geom_segment(
    x=17, xend=17, y=41, yend=39, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Right middle arrow
  ggplot2::geom_segment(
    x=83, xend=83, y=41, yend=39, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +  
  # Centre  lower arrow 
  ggplot2::geom_segment(
    x=50, xend=50, y=28, yend=26, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed")) +
  # Left lower arrow
  ggplot2::geom_segment(
    x=17, xend=17, y=28, yend=26, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type= "closed"))  
ggplot2::ggsave(p,width = 6, height = 8, units = 'in', dpi = 'retina', filename = "Flowchart.pdf", device='pdf', path = here::here("output"))
save(p, file = here::here("output","flowchart.RData"))