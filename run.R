# run.R
# renv::init()
# options(pkgType = "win.binary", install.packages.compile.from.source = "never")
# renv::snapshot(type = "all")

rm(list = ls()); gc(); Sys.setenv(RENV_CONFIG_GITIGNORE_READ = "FALSE"); options(repos = c(CRAN = "https://cran.rstudio.com/"), pkgType = "win.binary", install.packages.compile.from.source = "never")
# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", type = "win.binary"); try({ source("renv/activate.R") }, silent = TRUE); if(!requireNamespace("here", quietly=TRUE)) renv::install("here", type="win.binary", prompt=FALSE)
# if (file.exists("renv.lock")) { try({ renv::restore(prompt = FALSE) }, silent = TRUE); reqs <- gsub('.*"Package": "([^"]+)".*', "\\1", grep('"Package":', readLines("renv.lock", warn = FALSE), value = TRUE)); missing <- reqs[!reqs %in% dir(.libPaths()[1])]; if (length(missing) > 0) { message("Restoration failed for some packages. Force installing latest binaries..."); renv::install(missing, type = "win.binary", prompt = FALSE) } }

library(here)

# README --------------------------------------------------------------
# /data
#   ├ raw.csv      ← IH 추출 Chip 데이터 (칩 레벨)
#   └ ROOTID.csv   ← ROOTID - GROUP 매핑 (웨이퍼 레벨)
#
# /output
#   └ results.csv  ← DRB rev1.0 결과


# parameter -----------------------------------------------------------
# Only edit here!

# --- input files ---
FILE_NAME_RAW      <- "raw.csv"
FILE_NAME_META     <- "ROOTID.csv"

# --- key column names ---
ROOTID_COL         <- "ROOTID"
GROUP_COL          <- "GROUP"
PARTID_COL         <- "PARTID"
RADIUS_COL         <- "Radius"     # 실무: 좌/우 signed position (-150~150)

# --- group definition ---
# 원칙: REF=기준(Old), TARGET=변경(New)
GROUP_REF_LABEL     <- "A"
GROUP_TARGET_LABEL  <- "B"

# --- DRB rev1.0 parameters -------------------------------------------

# 1) Sigma_shift (k-sigma rule)
SIGMA_LEVEL         <- 1.0          # 예: 0.5 / 1.0 / 1.5

# 2) Wilcoxon (distribution shift)
WILCOX_ALPHA        <- 0.001         # 강한 신호만

# 3) KS (Radius-weighted KS)
# - ks_flag_radius_weighted() 사용
KS_ALPHA            <- 0.001         # 강한 신호만
KS_MIN_N            <- 30           # KS에 들어갈 최소 샘플 수(각 그룹)

# --- output ---
OUT_DIR             <- "output"
OUT_FILE_RESULTS    <- "results.csv"


# start ---------------------------------------------------------------
source(here::here("main.R"), encoding = "UTF-8")



