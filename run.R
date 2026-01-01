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

# --- key column names ---
ROOTID_COL          <- "ROOTID"
GROUP_COL           <- "GROUP"
PARTID_COL          <- "PARTID"
RADIUS_COL          <- "Radius"     # 실무: 좌/우 signed position (-150~150)

# --- group definition ---
# 원칙: REF=기준(Old), TARGET=변경(New)
GROUP_REF_LABEL     <- "A"
GROUP_TARGET_LABEL  <- "B"

# --- DRB rev1 parameters --------------------------------------------

# 1) direction 판정용 threshold (k-sigma)
# - sigma_score = (mean_target - mean_ref) / sd_ref
# - direction: Up/Down/Stable을 나누는 기준
SIGMA_LEVEL         <- 0.3          # 예: 0.3 / 0.5 / 1.0 / 1.5

# 2) spatial_drift (Sinkhorn OT) parameters
# Phase 1: clipping quantile (top 20%만 남기는 설정이면 0.8)
OT_Q                <- 0.80         # quantile threshold τ = Q(value, OT_Q)

# Phase 3: Sinkhorn regularization + iterations
OT_EPSILON          <- 0.10         # ε (regularization)
OT_MAX_ITER         <- 80           # 100 미만 권장 (속도)

# cost scaling (optional)
# - 좌표가 (0~300) 스케일이면 자동으로도 되는데,
#   데이터마다 스케일 다르면 여기서 조절 가능
OT_COST_SCALE       <- NULL         # NULL이면 내부에서 자동/기본 사용

# numeric stability
OT_TINY             <- 1e-12        # 0나눗셈 방지용

# --- output ---
OUT_DIR             <- "output"
OUT_FILE_RESULTS    <- "results.csv"

# start ---------------------------------------------------------------
source(here::here("main.R"), encoding = "UTF-8")










