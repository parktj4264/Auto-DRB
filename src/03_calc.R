# 03_calc.R
# src/03_calc.R ------------------------------------------------------------

cat("Calculating results (chip pooling)...\n")

# 0) 기본 체크 ---------------------------------------------------------------
if (!exists("dt")) stop("dt가 없습니다. 01_load.R 실행 여부를 확인하세요.")
if (!exists("MSR_COLS") || length(MSR_COLS) == 0) stop("MSR_COLS가 없습니다. 01_load.R의 MSR 식별을 확인하세요.")

need_cols <- c(ROOTID_COL, GROUP_COL, PARTID_COL)
miss_cols <- setdiff(need_cols, names(dt))
if (length(miss_cols) > 0) {
  stop(sprintf("dt에 필수 컬럼이 없습니다: %s", paste(miss_cols, collapse = ", ")))
}

# spatial_drift는 X/Y 필요
has_xy <- all(c(X_COL, Y_COL) %in% names(dt))
if (!has_xy) {
  cat("WARNING: dt에 X/Y 컬럼이 없습니다. spatial_drift는 NA로 저장됩니다.\n")
}

# REF/TARGET 라벨 체크
u_groups_dt <- sort(unique(as.character(dt[[GROUP_COL]])))
if (!(GROUP_REF_LABEL %in% u_groups_dt) || !(GROUP_TARGET_LABEL %in% u_groups_dt)) {
  cat("WARNING: dt의 GROUP 값에 REF/TARGET 라벨이 완전히 존재하지 않습니다.\n")
  cat(sprintf(" - dt GROUP levels: %s\n", paste(u_groups_dt, collapse = ", ")))
  cat(sprintf(" - REF/TARGET: %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
  cat("가능한 범위에서 계산을 시도합니다. (표본 부족 시 score는 NA 처리)\n")
}

# 1) 루프 준비 ---------------------------------------------------------------
n_msr <- length(MSR_COLS)
cat(sprintf(">> MSR count = %d\n", n_msr))



g_all <- as.character(dt[[GROUP_COL]])

idx_ref <- which(g_all == GROUP_REF_LABEL)
idx_target <- which(g_all == GROUP_TARGET_LABEL)

# OT spatial_drift 파라미터 (run.R에서 정의)
if (!exists("OT_SIGMA_THRESH")) OT_SIGMA_THRESH <- 3.0
if (!exists("OT_SMOOTH_SIGMA")) OT_SMOOTH_SIGMA <- 1.0
if (!exists("OT_MASK_GAIN")) OT_MASK_GAIN <- 2.0
if (!exists("OT_EPSILON")) OT_EPSILON <- 0.1
if (!exists("OT_MAX_ITER")) OT_MAX_ITER <- 80
if (!exists("OT_COST_SCALE")) OT_COST_SCALE <- NULL
if (!exists("OT_EMPTY_PENALTY")) OT_EMPTY_PENALTY <- 1.0
if (!exists("OT_TINY")) OT_TINY <- 1e-12

# 2) (핵심) GROUP 평균 맵을 미리 한 번에 만들어두기 ---------------------------

dtA_map <- NULL
dtB_map <- NULL

if (has_xy && length(idx_ref) > 0 && length(idx_target) > 0) {
  # cat(">> Precomputing GROUP mean maps (A/B) for all MSRs...\n")

  # 임시 좌표 컬럼 생성 (data.table inside-scope 안전하게)
  dt[, .x__ := suppressWarnings(as.numeric(get(X_COL)))]
  dt[, .y__ := suppressWarnings(as.numeric(get(Y_COL)))]

  ok_xy <- is.finite(dt$.x__) & is.finite(dt$.y__)

  # REF 맵
  dtA_map <- dt[
    ok_xy & get(GROUP_COL) == GROUP_REF_LABEL,
    lapply(.SD, mean, na.rm = TRUE),
    .SDcols = MSR_COLS,
    by = .(x = .x__, y = .y__)
  ]

  # TARGET 맵
  dtB_map <- dt[
    ok_xy & get(GROUP_COL) == GROUP_TARGET_LABEL,
    lapply(.SD, mean, na.rm = TRUE),
    .SDcols = MSR_COLS,
    by = .(x = .x__, y = .y__)
  ]

  # 정렬(선택)
  data.table::setorder(dtA_map, x, y)
  data.table::setorder(dtB_map, x, y)

  # 임시 컬럼 제거(깔끔 유지)
  dt[, c(".x__", ".y__") := NULL]

  # cat(sprintf("   -> dtA_map: %s rows, %s cols\n", format(nrow(dtA_map), big.mark=","), ncol(dtA_map)))
  # cat(sprintf("   -> dtB_map: %s rows, %s cols\n", format(nrow(dtB_map), big.mark=","), ncol(dtB_map)))
}


# 작은 헬퍼: 미리 계산된 맵에서 msr 컬럼만 v로 뽑기
get_map_from_precomputed <- function(map_dt, msr_col) {
  if (is.null(map_dt) || nrow(map_dt) == 0) {
    return(data.table::data.table(x = numeric(0), y = numeric(0), v = numeric(0)))
  }
  # .SDcols로 해당 MSR만 꺼내서 v로 이름 통일
  map_dt[, .(x, y, v = .SD[[1]]), .SDcols = msr_col]
}


# 3) MSR Loop (Parallel / Sequential) ----------------------------------------
if (!exists("N_CORES")) N_CORES <- 1

# --- Worker Function ---
calc_msr_single <- function(msr) {
  tryCatch(
    {
      # 1. Base Stats
      x_all <- dt[[msr]]

      sp <- split_ref_target(
        x = x_all,
        g = g_all,
        ref_label = GROUP_REF_LABEL,
        target_label = GROUP_TARGET_LABEL
      )

      x_ref <- sp$x_ref
      x_target <- sp$x_target

      # 요약통계
      summ <- calc_summary_ref_target(x_ref, x_target)

      # sigma_score
      sigma_score <- calc_sigma_score(
        mean_ref    = summ$mean_ref,
        sd_ref      = summ$sd_ref,
        mean_target = summ$mean_target
      )

      direction <- direction_from_sigma(
        sigma_score = sigma_score,
        sigma_level = SIGMA_LEVEL
      )

      cliffs_delta <- calc_cliffs_delta(x_ref, x_target)

      # spatial_drift
      spatial_drift <- NA_real_
      if (!is.null(dtA_map) && !is.null(dtB_map) && nrow(dtA_map) > 0 && nrow(dtB_map) > 0) {
        map_ref <- get_map_from_precomputed(dtA_map, msr)
        map_target <- get_map_from_precomputed(dtB_map, msr)

        spatial_drift <- calc_spatial_drift_sinkhorn(
          map_ref       = map_ref,
          map_target    = map_target,
          sigma_thresh  = OT_SIGMA_THRESH,
          smooth_sigma  = OT_SMOOTH_SIGMA,
          mask_gain     = OT_MASK_GAIN,
          epsilon       = OT_EPSILON,
          max_iter      = OT_MAX_ITER,
          cost_scale    = OT_COST_SCALE,
          empty_penalty = OT_EMPTY_PENALTY,
          tiny          = OT_TINY
        )
      }

      # Return Row
      make_result_row(
        msr_name      = msr,
        direction     = direction,
        sigma_score   = sigma_score,
        cliffs_delta  = cliffs_delta,
        spatial_drift = spatial_drift,
        mean_ref      = summ$mean_ref,
        sd_ref        = summ$sd_ref,
        mean_target   = summ$mean_target,
        sd_target     = summ$sd_target
      )
    },
    error = function(e) {
      # Error Handling: Return NA row with warning
      # cat(sprintf("\n[ERROR] MSR: %s -> %s\n", msr, e$message)) # Parallel에서 cat은 잘 안보일 수 있음

      make_result_row(
        msr_name      = msr,
        direction     = "Error",
        sigma_score   = NA_real_,
        cliffs_delta  = NA_real_,
        spatial_drift = NA_real_,
        mean_ref      = NA_real_,
        sd_ref        = NA_real_,
        mean_target   = NA_real_,
        sd_target     = NA_real_
      )
    }
  )
}

# --- Execution ---
if (N_CORES > 1) {
  cat(sprintf(">> Setting up Parallel Mode (Cores: %d)...\n", N_CORES))
  future::plan(future::multisession, workers = N_CORES)
} else {
  cat(">> Running in Sequential Mode\n")
  future::plan(future::sequential)
}

cat(">> Processing MSRs...\n")

# Use progressr for progress bar
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

res_list <- progressr::with_progress({
  p <- progressr::progressor(along = MSR_COLS)

  future.apply::future_lapply(
    MSR_COLS,
    function(msr) {
      res <- calc_msr_single(msr)
      p() # update progress
      res
    },
    future.seed = TRUE
  )
})

# Reset plan
future::plan(future::sequential)

# 4) 결과 합치기 --------------------------------------------------------------
results_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)

# 5) 정렬 (권장: abs(sigma_score) 내림차순)
if ("sigma_score" %in% names(results_dt)) {
  ord <- order(is.na(results_dt$sigma_score), -abs(results_dt$sigma_score))
  results_dt <- results_dt[ord]
}

cat("Done (Calc)\n")
cat(sprintf(
  "Results: %s rows, %s cols\n",
  format(nrow(results_dt), big.mark = ","), ncol(results_dt)
))

if (all(c("direction") %in% names(results_dt))) {
  tab_dir <- table(results_dt$direction, useNA = "ifany")
  cat("Direction summary:\n")
  print(tab_dir)
}
