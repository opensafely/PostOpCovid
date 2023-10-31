library(foreach)
library(data.table)
library(magrittr)
library(survival)

ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis","Utils.R"))

## Load data files
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output","cohort_long.feather")))

## Convert procedures to individual categorical
dt.tv[, procedure := ifelse(Abdominal, "Abdominal",
                            ifelse(Obstetrics, "Obstetrics",
                                   ifelse(Orthopaedic, "Orthopaedic",
                                          ifelse(CardioThoracicVascular, "CardioThoracicVascular", NA))))]

# Set Abdominal as reference
dt.tv[, procedure := as.factor(procedure)]
dt.tv[, procedure := relevel(procedure, ref = "Abdominal")]


## Define covariates
covariates <- c('procedure','age.cat','sex','bmi.cat','imd5','postcovid','wave','Major.op',
                'vaccination.status.factor','region','Current.Cancer','Emergency','LOS.bin','Charl12','recentCOVID','previousCOVID')

data.table::setkey(dt.tv,patient_id,tstart,tstop)
# set region as factor 
dt.tv[, region := as.factor(region)]


# Functions for crude HR and Cumulative Incidence for VTE 
crude.HR.vte <- function(x) {
  vapply(1:length(coef(x)), FUN.VALUE = '', function(i) paste0(round(exp(coef(x)[i]),2), ' (', paste(round(exp(confint(x)[i,]),2), collapse = ','),')'))
}


#Function for adjusted HR function
adjusted.HR.vte <- function(covariate_name) {
  fully_adjusted_model <- coxph(formula = as.formula(paste0('Surv(start,end,event.VTE==1) ~ ', paste(covariates, collapse = " + "))),
                                id = patient_id,
                                data = dt.tv[start >= 0 & any.op.VTE == T])
  
  coef_name <- paste0(covariate_name, collapse = ", ")
  
  if (all(coef_name %in% names(coef(fully_adjusted_model)))) {
    coef_val <- coef(fully_adjusted_model)[coef_name]
    conf_int_vals <- as.matrix(confint(fully_adjusted_model)[coef_name, ])
    return(paste0(round(exp(coef_val), 2), " (", round(exp(conf_int_vals[1]), 2), ",", round(exp(conf_int_vals[2]), 2), ")"))
  } else {
    return(rep(NA, length(coef_name)))
  }
}

# for(i in 1:length(covariates)){
#   print (covariates[i])
#   print(dt.tv[start >= 0 ,.N, event.VTE][])
#   dt.tv[start >= 0 & (tstop - tstart) <= 30 ,.N, event.VTE][]
#   dt.tv[start >= 0 & (tstop - tstart) <= 30 & any.op.VTE == T ,.N, event.VTE][]
#   dt.tv[start >= 0 & (tstop - tstart) <= 30 & any.op.VTE == T & !is.na(get (covariates[i])),.N, event.VTE][]
  
# }
## Crude and fully adjusted HR calculations
crude.vte.cov <- data.table::rbindlist(
  lapply(1:length(covariates), function(i) {
    
    
    # Define covariates, level, cumulative incidences, crude and adjusted HR with 95% CI- NB- any.op.VTE==T replaced with any.op==T 
    covariate_names <- rep(covariates[i], length(levels(dt.tv[start >=0 & any.op == T, as.factor(get(covariates[i]))])))
    
    covariate_levels <- levels(dt.tv[start >=0 & any.op == T, as.factor(get(covariates[i]))])
    
    cuminc_values <- cuminc.km.vte(covariates[i], niter = 2)[,2:4]
    
    crude_hr_values <- c('Ref', coxph(formula = as.formula(paste0('Surv(start,end,event.VTE==1) ~ ', covariates[i])),
                                      id = patient_id,
                                      data = dt.tv[start >=0 & any.op == T]) %>% crude.HR.vte())
    
    adjusted_hr_values <- sapply(covariate_levels, function(lev) {
      adjusted.HR.vte(paste0(covariates[i], lev))
    })
    
    cbind(covariate_names, covariate_levels, cuminc_values, crude_hr_values, adjusted_hr_values)
  })
)



# Name the columns of the output table
names(crude.vte.cov) <- c("Characteristic", "Level", "Number at risk", "Number of post operative VTE events within 30 days", 
                          "30 day Cumulative VTE Risk adjusted for censoring", "Crude Hazard ratio 30 day VTE (95% CI)", "Adjusted Hazard ratio 30 day VTE (95% CI)")

# Save the results to csv and RData
data.table::fwrite(crude.vte.cov, file = here::here("output","postopvte_fullyadjusted30dayanalysis.csv"))
save(crude.vte.cov, file = here::here("output","postopvte_fullyadjusted30dayanalysis.RData"))

##################################




