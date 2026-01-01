# 03_calc.R
# src/03_calc.R ------------------------------------------------------------

cat("Calculating results (chip pooling)...\n")

# 0) 기본 체크 ---------------------------------------------------------------
if (!exists("dt")) stop("dt가 없습니다. 01_load.R 실행 여부를 확인하세요.")
if (!exists("MSR_COLS") || length(MSR_COLS) == 0) stop("MSR_COLS가 없습니다. 01_load.R의 MSR 식별을 확인하세요.")

# 필수 컬럼 존재 체크 (경고 수준으로 최대한 유연하게)
need_cols <- c(ROOTID_COL, GROUP_COL, RADIUS_COL, PARTID_COL)
miss_cols <- setdiff(need_cols, names(dt))
if (length(miss_cols) > 0) {
  stop(sprintf("dt에 필수 컬럼이 없습니다: %s", paste(miss_cols, collapse=", ")))
}

# OT spatial_drift는 X/Y가 필요 (없으면 spatial_drift만 NA로 처리)
has_xy <- all(c("X","Y") %in% names(dt))
if (!has_xy) {
  cat("WARNING: dt에 X/Y 컬럼이 없습니다. spatial_drift는 NA로 저장됩니다.\n")
}

# REF/TARGET 라벨 체크 (유연: 없으면 경고 후 가능한 만큼)
u_groups_dt <- sort(unique(as.character(dt[[GROUP_COL]])))
if (!(GROUP_REF_LABEL %in% u_groups_dt) || !(GROUP_TARGET_LABEL %in% u_groups_dt)) {
  cat("WARNING: dt의 GROUP 값에 REF/TARGET 라벨이 완전히 존재하지 않습니다.\n")
  cat(sprintf(" - dt GROUP levels: %s\n", paste(u_groups_dt, collapse=", ")))
  cat(sprintf(" - REF/TARGET: %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
  cat("가능한 범위에서 계산을 시도합니다. (표본 부족 시 score는 NA 처리)\n")
}

# 1) 루프 준비 ---------------------------------------------------------------
n_msr <- length(MSR_COLS)
cat(sprintf(">> MSR count = %d\n", n_msr))

# 진행바
pb <- utils::txtProgressBar(min = 0, max = n_msr, style = 3)

# 결과 저장 (리스트로 모아 rbindlist)
res_list <- vector("list", n_msr)

# 미리 뽑아두는 벡터 (속도)
g_all <- as.character(dt[[GROUP_COL]])
r_all <- dt[[RADIUS_COL]]  # 여기선 아직 사용 안 하더라도 (다른 지표 확장 대비) 유지

# 그룹 인덱스는 루프 밖에서 한 번만 (속도 + 일관성)
idx_ref    <- which(g_all == GROUP_REF_LABEL)
idx_target <- which(g_all == GROUP_TARGET_LABEL)

# OT spatial_drift 파라미터 (run.R에서 정의)
if (!exists("OT_Q")) OT_Q <- 0.8
if (!exists("OT_EPSILON")) OT_EPSILON <- 0.1
if (!exists("OT_MAX_ITER")) OT_MAX_ITER <- 80
if (!exists("OT_COST_SCALE")) OT_COST_SCALE <- NULL
if (!exists("OT_TINY")) OT_TINY <- 1e-12

# 2) MSR loop ----------------------------------------------------------------
for (i in seq_len(n_msr)) {
  
  msr <- MSR_COLS[i]
  
  # MSR 벡터
  x_all <- dt[[msr]]
  
  # REF/TARGET split (칩 풀링)
  sp <- split_ref_target(
    x = x_all,
    g = g_all,
    ref_label = GROUP_REF_LABEL,
    target_label = GROUP_TARGET_LABEL
  )
  
  x_ref    <- sp$x_ref
  x_target <- sp$x_target
  
  # 요약통계
  summ <- calc_summary_ref_target(x_ref, x_target)
  
  # sigma_score (Glass's delta)
  sigma_score <- calc_sigma_score(
    mean_ref    = summ$mean_ref,
    sd_ref      = summ$sd_ref,
    mean_target = summ$mean_target
  )
  
  # direction (Up/Down/Stable) : threshold = SIGMA_LEVEL
  direction <- direction_from_sigma(
    sigma_score  = sigma_score,
    sigma_level  = SIGMA_LEVEL
  )
  
  # cliffs_delta (distribution dominance)
  cliffs_delta <- calc_cliffs_delta(
    x_ref    = x_ref,
    x_target = x_target
  )
  
  # spatial_drift (Sinkhorn OT on GROUP-averaged wafer maps)
  spatial_drift <- NA_real_
  if (has_xy && length(idx_ref) > 0 && length(idx_target) > 0) {
    
    # 그룹별 좌표-평균 맵 생성 (칩 좌표별 평균)
    map_ref <- make_group_mean_map(
      dt          = dt,
      msr_col     = msr,
      group_col   = GROUP_COL,
      group_label = GROUP_REF_LABEL,
      x_col       = "X",
      y_col       = "Y"
    )
    
    map_target <- make_group_mean_map(
      dt          = dt,
      msr_col     = msr,
      group_col   = GROUP_COL,
      group_label = GROUP_TARGET_LABEL,
      x_col       = "X",
      y_col       = "Y"
    )
    
    spatial_drift <- calc_spatial_drift_sinkhorn(
      map_ref     = map_ref,
      map_target  = map_target,
      q_clip      = OT_Q,
      epsilon     = OT_EPSILON,
      max_iter    = OT_MAX_ITER,
      cost_scale  = OT_COST_SCALE,
      tiny        = OT_TINY
    )
  }
  
  # 결과 row 생성 (Final schema)
  res_list[[i]] <- make_result_row(
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
  
  # progress update
  utils::setTxtProgressBar(pb, i)
}

close(pb)

# 3) 결과 합치기 --------------------------------------------------------------
results_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)

# 4) 정렬 (권장: abs(sigma_score) 내림차순) -----------------------------------
# - NA는 뒤로
if ("sigma_score" %in% names(results_dt)) {
  ord <- order(is.na(results_dt$sigma_score), -abs(results_dt$sigma_score))
  results_dt <- results_dt[ord]
}

cat("Done (Calc)\n")
cat(sprintf("Results: %s rows, %s cols\n",
            format(nrow(results_dt), big.mark=","), ncol(results_dt)))

# quick summary
if (all(c("direction") %in% names(results_dt))) {
  tab_dir <- table(results_dt$direction, useNA = "ifany")
  cat("Direction summary:\n")
  print(tab_dir)
}
