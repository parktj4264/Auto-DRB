# run.R
# renv::init()
# options(pkgType = "win.binary", install.packages.compile.from.source = "never")
# renv::snapshot(type = "all")

rm(list = ls()); gc(); Sys.setenv(RENV_CONFIG_GITIGNORE_READ = "FALSE"); options(repos = c(CRAN = "https://cran.rstudio.com/"), pkgType = "win.binary", install.packages.compile.from.source = "never")
# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", type = "win.binary"); try({ source("renv/activate.R") }, silent = TRUE); if(!requireNamespace("here", quietly=TRUE)) renv::install("here", type="win.binary", prompt=FALSE)
# if (file.exists("renv.lock")) { try({ renv::restore(prompt = FALSE) }, silent = TRUE); reqs <- gsub('.*"Package": "([^"]+)".*', "\\1", grep('"Package":', readLines("renv.lock", warn = FALSE), value = TRUE)); missing <- reqs[!reqs %in% dir(.libPaths()[1])]; if (length(missing) > 0) { message("Restoration failed for some packages. Force installing latest binaries..."); renv::install(missing, type = "win.binary", prompt = FALSE) } }

library(here)


# READMD ------------------------------------------------------------------
# /data
#   └ raw.csv      ← IH 추출 Chip 데이터 넣기
# /output
#   ├ output1.csv
#   └ output2.csv


# parameter ---------------------------------------------------------------
# Only edit here!

FILE_NAME_RAW  <- "raw.csv"
FILE_NAME_META <- "ROOTID.csv"

SIGMA_LEVEL    <- 1.0            # Outlier 기준 (1 sigma)
WILCOX_MU      <- 0
WAFER_RADIUS   <- 150

ADD_GROUPS     <- c("공정", "CORE")
GROUP_COLS     <- c("ROOTID", "GROUP", ADD_GROUPS)


# start -------------------------------------------------------------------
source(here::here("main.R"), encoding = "UTF-8")
