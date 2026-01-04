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
#   └ results.csv  ← DRB rev1 결과


# parameter -----------------------------------------------------------
# Only edit here!

# --- input files ---
FILE_NAME_RAW       <- "raw.csv"
FILE_NAME_META      <- "ROOTID.csv"

# --- group definition ---
# 원칙: REF=기준(Old), TARGET=변경(New)
GROUP_REF_LABEL     <- "A"
GROUP_TARGET_LABEL  <- "B"


ROOTID_COL <- "ROOTID"
PARTID_COL <- "PARTID"
GROUP_COL  <- "GROUP"
X_COL <- "X"
Y_COL <- "Y"

# --- DRB rev1 parameters --------------------------------------------

# 1) direction 판정용 threshold (k-sigma)
# - sigma_score = (mean_target - mean_ref) / sd_ref
# - direction: Up/Down/Stable을 나누는 기준
SIGMA_LEVEL         <- 0.3          # 예: 0.3 / 0.5 / 1.0 / 1.5

# 2) spatial_drift (Robust preprocess -> Smooth -> Sinkhorn OT)
# Phase 1: Z-filter (abs(z) > thresh 만 에너지로 인정)
OT_SIGMA_THRESH     <- 3.0
# Phase 2: smoothing sigma (blob화)
OT_SMOOTH_SIGMA     <- 1.0
# Phase 3: Sinkhorn regularization + iterations
OT_EPSILON          <- 0.10
OT_MAX_ITER         <- 80
# cost scaling (optional)
OT_COST_SCALE       <- NULL  # NULL이면 내부에서 자동
# empty case handling
OT_EMPTY_PENALTY    <- 1.0   # 한쪽만 신호 있으면 penalty
# numeric stability
OT_TINY             <- 1e-12


# --- output ---
OUT_DIR             <- "output"
OUT_FILE_RESULTS    <- "results.csv"

# start ---------------------------------------------------------------
source(here::here("main.R"), encoding = "UTF-8")










