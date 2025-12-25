# main.R

# library load and setting ------------------------------------------------
library(here)
source(here::here("src", "00_libraries.R"), encoding = "UTF-8")

cat("==============================================\n")

# 시간 기록 시작
start_time <- Sys.time()


# function & script load --------------------------------------------------
source(here::here("src", "01_load.R"), encoding = "UTF-8")
source(here::here("src", "02_funcs.R"), encoding = "UTF-8")
source(here::here("src", "03_calc.R"), encoding = "UTF-8")
source(here::here("src", "04_save.R"), encoding = "UTF-8")



# 실행 완료 후 시간 출력 ---------------------------------------------------
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
mins <- floor(elapsed / 60)
secs <- round(elapsed %% 60)

cat(sprintf("\nElapsed time: %d min %d sec\n", mins, secs))

cat("==============================================\n")