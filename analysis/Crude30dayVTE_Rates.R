library(data.table)
ncores <- parallel::detectCores(logical = T) - 1
data.table::setDTthreads(ncores)

source(here::here("analysis", "Utils.R"))

# Load Data
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output", "cohort_long.feather")))
procedures <- c('Abdominal', 'Obstetrics', 'Orthopaedic', 'CardioThoracicVascular')
covariates <- c(procedures, 'age.cat', 'sex', 'bmi.cat', 'imd5', 'wave', 'Major.op', 'Intermediate.or.Minor.op', 'vaccination.status.factor', 'region', 'Current.Cancer', 'Emergency', 'LOS.bin', 'Charl12', 'recentCOVID', 'previousCOVID', 'postcovid')


data.table::setkey(dt.tv, patient_id, tstart, tstop)

summary(dt.tv)

# set region as factor
dt.tv[, region := as.factor(region)]

# Calculate Crude VTE Rates
dt.tv[, vte := max(as.numeric(event.VTE == 1), na.rm = TRUE), keyby = patient_id]

# Loop for Major/Intermediate/Minor and Emergency/Non-Emergency
for (op_type in c('Major.op', 'Intermediate.or.Minor.op')) {
  for (Emerg in c(F, T)) {
    Emerg_label <- ifelse(Emerg, "Emergency", "Elective")
    print(paste0(op_type,Emerg))
    crude.vte.rates <- cbind(
      rbindlist(
        lapply(covariates, function(x) {
          dt.tv[Emergency == Emerg & get(op_type) == T & (postcovid.VTE.cohort) & !is.na(get(x)) & start >= 0 &(tstop - tstart) <=30, head(.SD, 1), keyby = .(patient_id, end.fu)][
            , .N, keyby = c(x, 'vte')
          ][, .(Total = sum(N), Events = sum(N[vte == 1])), keyby = x][, .(covariate = x, levels = get(x), Total, Events, Crude.Rate = 100 * Events / Total)]
        })
      )
    )
    
    data.table::fwrite(crude.vte.rates, file = here::here("output", paste0(op_type, "_", Emerg_label, "30_day_Crude_VTE_Rates.csv")))
    #   save(crude.vte.rates, file = here::here("output", paste0(op_type, "_", Emerg_label, ".RData")))
  }
}