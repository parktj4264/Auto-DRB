# 03_calc.R

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

# REF/TARGET 라벨이 dt에 실제로 존재하는지 재확인 (유연: 없으면 경고 후 가능한 만큼)
u_groups_dt <- sort(unique(as.character(dt[[GROUP_COL]])))
if (!(GROUP_REF_LABEL %in% u_groups_dt) || !(GROUP_TARGET_LABEL %in% u_groups_dt)) {
  cat("WARNING: dt의 GROUP 값에 REF/TARGET 라벨이 완전히 존재하지 않습니다.\n")
  cat(sprintf(" - dt GROUP levels: %s\n", paste(u_groups_dt, collapse=", ")))
  cat(sprintf(" - REF/TARGET: %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
  cat("가능한 범위에서 계산을 시도합니다. (표본 부족 시 flag는 FALSE 처리)\n")
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

# 2) MSR loop ----------------------------------------------------------------
for (i in seq_len(n_msr)) {
  
  msr <- MSR_COLS[i]
  
  # MSR 벡터
  x_all <- dt[[msr]]
  
  # REF/TARGET split
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
  
  # Sigma flags
  sig <- sigma_shift_flags(
    mean_ref    = summ$mean_ref,
    sd_ref      = summ$sd_ref,
    mean_target = summ$mean_target,
    sd_target   = summ$sd_target,
    sigma_level = SIGMA_LEVEL
  )
  
  # Wilcox flag
  w <- wilcox_flag(
    x_ref    = x_ref,
    x_target = x_target,
    alpha    = WILCOX_ALPHA
  )
  
  # KS flag (Radius-weighted KS)
  # radius는 그룹별로 같은 row index로 split해야 하므로, 인덱스 기반으로 구성
  idx_ref    <- which(g_all == GROUP_REF_LABEL)
  idx_target <- which(g_all == GROUP_TARGET_LABEL)
  
  ks <- ks_flag_radius_weighted(
    radius_ref    = r_all[idx_ref],
    x_ref         = x_all[idx_ref],
    radius_target = r_all[idx_target],
    x_target      = x_all[idx_target],
    alpha         = KS_ALPHA,
    min_n         = KS_MIN_N
  )
  
  # 결과 row 생성
  res_list[[i]] <- make_result_row(
    msr_name     = msr,
    sigma_up     = sig$sigma_up,
    sigma_down   = sig$sigma_down,
    wilcox_flag  = w$flag,
    ks_flag      = ks$flag,
    mean_ref     = summ$mean_ref,
    sd_ref       = summ$sd_ref,
    mean_target  = summ$mean_target,
    sd_target    = summ$sd_target,
    mean_diff    = summ$mean_diff,
    median_diff  = summ$median_diff
  )
  
  # progress update
  utils::setTxtProgressBar(pb, i)
}

close(pb)

# 3) 결과 합치기 --------------------------------------------------------------
results_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)

cat("Done (Calc)\n")
cat(sprintf("Results: %s rows, %s cols\n",
            format(nrow(results_dt), big.mark=","), ncol(results_dt)))

# quick summary
cat(sprintf("Flags count: sigma_up=%d, sigma_down=%d, wilcox=%d, ks=%d\n",
            sum(results_dt$sigma_up, na.rm=TRUE),
            sum(results_dt$sigma_down, na.rm=TRUE),
            sum(results_dt$wilcox_flag, na.rm=TRUE),
            sum(results_dt$ks_flag, na.rm=TRUE)
))
