# Load libraries

library(data.table)  
library(survival)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load data
dt.tv <- data.table::setDT(arrow::read_feather(here::here("output", "cohort_long.feather")))

#############################################

# Define procedures and waves
procedures <- c("Abdominal", "Obstetrics", "Orthopaedic", "CardioThoracicVascular")
waves <- c("Wave_1", "Wave_2", "Wave_3", "Wave_4")

# Define a function to calculate outcomes for each procedure and wave

cuminc.km.vte_wave <- function(proc, wv) {
  safelog <- function(x) {
    x[x < 1e-200] <- 1e-200
    log(x)
  }
  data.table::setkey(dt.tv, patient_id, tstart, tstop)
  
  # Filter the dataset for the current procedure and wave
  filtered_data <- dt.tv[get(proc) == 1 & wave == wv & start >= 0 & (tstop - tstart) <= 30 & any.op.VTE == TRUE]
  
  # Check if there are observations in the filtered data
  if (nrow(filtered_data) > 0) {
    
    # Perform survival analysis to calculate outcomes
    surv_fit <- survival::survfit(survival::Surv(start, end, event.VTE == 1) ~ 1, data = filtered_data, id = patient_id)
    n_risk <- summary(surv_fit, times = 0)$n.risk
    n_event <- summary(surv_fit, times = 30)$n.event
    surv_30_day <- 1 - summary(surv_fit, times = 30)$surv
    cuminc_30_day_incidence <- 100 * round(surv_30_day, digits = 4)
  } else {
    # If no data, set outcomes to 0
    n_risk <- 0
    n_event <- 0
    cuminc_30_day_incidence <- 0
  }
  
  # Return the calculated outcomes
  return(data.table(Procedure = proc, Wave = wv, num_at_risk = n_risk, num_events = n_event, cuminc_30_day = cuminc_30_day_incidence))
}

# Create an empty results table to store the outcomes
procedure_across_waves <- data.table(Procedure = character(),
                           Wave = character(),
                           num_at_risk = numeric(),
                           num_events = numeric(),
                           cuminc_30_day = numeric())
  
# Calculate outcomes for each combination of procedure and wave
for (proc in procedures) {
  for (wv in waves) {
    procedure_across_waves <- rbind(procedure_across_waves, cuminc.km.vte_wave(proc, wv))
  }
}

# Check output in console
print(procedure_across_waves)

# Name the columns of the results table
names(procedure_across_waves) <- c("Procedure", "Wave", "Number at risk", "Number of post operative VTE events within 30 days", 
                    "30 day Cumulative VTE")

# Save the results to CSV and RData
data.table::fwrite(procedure_across_waves, file = here::here("output", "table_postop_VTEincidence_by_procedure_wave.csv"))
save(procedure_across_waves, file = here::here("output", "table_postop_VTEincidence_by_procedure_wave.RData"))

### Plot 30 day incidence across waves

VTE_30dayincidence_procedure_waves_plot <- ggplot2::ggplot(procedure_across_waves, aes(x = Wave, 
                                                              y = `30 day Cumulative VTE`, 
                                                              group = Procedure,
                                                              color = Procedure)) +
  ggplot2::geom_line() +
  ggplot2::labs(title = "30-Day Cumulative VTE Incidence Across Procedures and Waves",
                x = "Wave",
                y = "30-Day Cumulative VTE Incidence (%)") +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"))  

# Save the plot as an image file
ggplot2::ggsave(plot = VTE_waves_plot, 
                filename = here::here("output", "VTE_waves_plot.png"), dpi = 300, width = 7, height = 5, units = "in", device = "png")


