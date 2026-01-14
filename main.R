# main.R -------------------------------------------------------------------

# library load and setting --------------------------------------------
library(here)
source(here::here("src", "00_libraries.R"), encoding = "UTF-8")

# 0) Logging Setup ----------------------------------------------------
if (!dir.exists("log")) dir.create("log")
log_file <- sprintf("log/log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))
sink(log_file, split = TRUE) # Console + File

cat("==============================================\n")
cat(sprintf("Log File: %s\n", log_file))

# time start -----------------------------------------------------------
start_time <- Sys.time()

# function & script load ----------------------------------------------
source(here::here("src", "01_load.R"), encoding = "UTF-8")
source(here::here("src", "02_funcs.R"), encoding = "UTF-8")
source(here::here("src", "03_calc.R"), encoding = "UTF-8")
source(here::here("src", "04_save.R"), encoding = "UTF-8")

# time end -------------------------------------------------------------
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
mins <- floor(elapsed / 60)
secs <- round(elapsed %% 60)

cat(sprintf("\nElapsed time: %d min %d sec\n", mins, secs))
cat("DRB rev1 done\n")
cat("==============================================\n")

# Close Logging
sink()
