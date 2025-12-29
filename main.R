# main.R -------------------------------------------------------------------

# library load and setting --------------------------------------------
library(here)
source(here::here("src", "00_libraries.R"), encoding = "UTF-8")

cat("==============================================\n")

# Basic run info -------------------------------------------------------
cat("----------------------------------------------\n")
cat("DRB rev1 start\n")
cat(sprintf("- raw file            : %s\n", FILE_NAME_RAW))
cat(sprintf("- meta file           : %s\n", FILE_NAME_META))
cat(sprintf("- group (REF/TARGET)  : %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
cat(sprintf("- sigma level         : %s\n", SIGMA_LEVEL))
cat(sprintf("- ws bins (K)         : %s\n", WS_N_BINS))
cat(sprintf("- ws bin method       : %s\n", WS_BIN_METHOD))
cat(sprintf("- ws min n / bin      : %s\n", WS_MIN_N_PER_BIN))
cat(sprintf("- output              : %s/%s\n", OUT_DIR, OUT_FILE_RESULTS))
cat("----------------------------------------------\n")

# time start -----------------------------------------------------------
start_time <- Sys.time()

# function & script load ----------------------------------------------
source(here::here("src", "01_load.R"), encoding = "UTF-8")
source(here::here("src", "02_funcs.R"), encoding = "UTF-8")
source(here::here("src", "03_calc.R"), encoding = "UTF-8")
source(here::here("src", "04_save.R"), encoding = "UTF-8")

# time end -------------------------------------------------------------
end_time <- Sys.time()
elapsed  <- as.numeric(difftime(end_time, start_time, units = "secs"))
mins     <- floor(elapsed / 60)
secs     <- round(elapsed %% 60)

cat(sprintf("\nElapsed time: %d min %d sec\n", mins, secs))
cat("DRB rev1 done\n")
cat("==============================================\n")
