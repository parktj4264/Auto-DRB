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
r_all <- dt[[RADIUS_COL]]

# 그룹 인덱스는 루프 밖에서 한 번만 (속도 + 일관성)
idx_ref    <- which(g_all == GROUP_REF_LABEL)
idx_target <- which(g_all == GROUP_TARGET_LABEL)

# ws_spatial 파라미터 (run.R에서 정의)
if (!exists("WS_N_BINS")) WS_N_BINS <- 30
if (!exists("WS_BIN_METHOD")) WS_BIN_METHOD <- "equal_width"
if (!exists("WS_MIN_N_PER_BIN")) WS_MIN_N_PER_BIN <- 10

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
  
  # ws_spatial (radius profile L1 deviation)
  # - radius는 같은 row index로 그룹별 분리해야 하므로 idx_ref/idx_target 사용
  ws_spatial <- calc_ws_spatial(
    radius_ref    = r_all[idx_ref],
    x_ref         = x_all[idx_ref],
    radius_target = r_all[idx_target],
    x_target      = x_all[idx_target],
    K             = WS_N_BINS,
    method        = WS_BIN_METHOD,
    min_n_per_bin = WS_MIN_N_PER_BIN
  )
  
  # 결과 row 생성 (Final schema)
  res_list[[i]] <- make_result_row(
    msr_name      = msr,
    direction     = direction,
    sigma_score   = sigma_score,
    cliffs_delta  = cliffs_delta,
    ws_spatial    = ws_spatial,
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
