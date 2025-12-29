# 02_funcs.R
# src/02_funcs.R ------------------------------------------------------------

# 목적:
# - 03_calc.R에서 MSR 루프를 돌릴 때 사용할 함수 모음
# - REF/TARGET 비교 기준은 run.R의 GROUP_REF_LABEL / GROUP_TARGET_LABEL
# - 결과는 results.csv (snake_case, ref/target 표기) 기준
#
# DRB rev1 (Final):
# - sigma_score      : Glass's delta = (mean_target - mean_ref) / sd_ref
# - cliffs_delta     : distribution dominance (Wilcoxon W -> U -> delta)
# - ws_spatial       : Radius profile L1 deviation (bin mean profile gap)
# - direction        : Up/Down/Stable (threshold = SIGMA_LEVEL)

# -------------------------------------------------------------------------
# 0) Small helpers
# -------------------------------------------------------------------------

safe_mean <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

# 그룹 split (벡터 기반)
split_ref_target <- function(x, g, ref_label, target_label) {
  idx_ref    <- which(g == ref_label)
  idx_target <- which(g == target_label)
  
  list(
    x_ref    = x[idx_ref],
    x_target = x[idx_target],
    n_ref    = length(idx_ref),
    n_target = length(idx_target)
  )
}

# -------------------------------------------------------------------------
# 1) Summary stats (chip pooling)
# -------------------------------------------------------------------------

calc_summary_ref_target <- function(x_ref, x_target) {
  mean_ref    <- safe_mean(x_ref)
  sd_ref      <- safe_sd(x_ref)
  mean_target <- safe_mean(x_target)
  sd_target   <- safe_sd(x_target)
  
  mean_diff   <- mean_target - mean_ref
  median_diff <- safe_median(x_target) - safe_median(x_ref)
  
  list(
    mean_ref    = mean_ref,
    sd_ref      = sd_ref,
    mean_target = mean_target,
    sd_target   = sd_target,
    mean_diff   = mean_diff,
    median_diff = median_diff
  )
}

# -------------------------------------------------------------------------
# 2) sigma_score (Glass's delta)
#    score = (mean_target - mean_ref) / sd_ref
# -------------------------------------------------------------------------

calc_sigma_score <- function(mean_ref, sd_ref, mean_target) {
  if (!is.finite(mean_ref) || !is.finite(mean_target) || !is.finite(sd_ref)) return(NA_real_)
  if (sd_ref <= 0) return(NA_real_)
  (mean_target - mean_ref) / sd_ref
}

direction_from_sigma <- function(sigma_score, sigma_level) {
  # sigma_level = SIGMA_LEVEL (예: 0.5 / 1.0 / 1.5)
  if (!is.finite(sigma_score) || !is.finite(sigma_level)) return("Stable")
  if (sigma_score >= sigma_level) return("Up")
  if (sigma_score <= -sigma_level) return("Down")
  "Stable"
}

# -------------------------------------------------------------------------
# 3) cliffs_delta (Stochastic Dominance)
#    - Wilcoxon rank-sum의 W 통계량으로부터 U를 복원
#    - delta = 2U/(n_t*n_r) - 1  (Target > Ref 우세면 +)
#
#    주의:
#    - R의 wilcox.test(x, y)$statistic 은 "첫 번째 샘플 x"의 rank sum(W)임.
#    - 따라서 Target 우세를 +로 만들려면 wilcox.test(x_target, x_ref)로 호출.
# -------------------------------------------------------------------------

calc_cliffs_delta <- function(x_ref, x_target) {
  xr <- x_ref[is.finite(x_ref)]
  xt <- x_target[is.finite(x_target)]
  
  n_r <- length(xr)
  n_t <- length(xt)
  
  if (n_r < 1 || n_t < 1) return(NA_real_)
  
  W <- tryCatch({
    as.numeric(stats::wilcox.test(xt, xr, alternative = "two.sided", exact = FALSE)$statistic)
  }, error = function(e) NA_real_)
  
  if (!is.finite(W)) return(NA_real_)
  
  # U = W - n_t(n_t+1)/2
  U <- W - (n_t * (n_t + 1)) / 2
  
  denom <- n_t * n_r
  if (denom <= 0) return(NA_real_)
  
  delta <- (2 * U / denom) - 1
  as.numeric(delta)
}

# -------------------------------------------------------------------------
# 4) ws_spatial (Radius profile L1 deviation)
#    - Radius를 K bins로 나눔
#    - 각 bin에서 local mean profile 생성
#    - score = mean( abs(profile_ref - profile_target) ) over valid bins
#
#    params:
#      K             : WS_N_BINS
#      method        : WS_BIN_METHOD ("equal_width" or "quantile")
#      min_n_per_bin : WS_MIN_N_PER_BIN (각 그룹 bin 최소 표본)
# -------------------------------------------------------------------------

make_radius_bins <- function(r_all, K, method = "equal_width") {
  r_all <- r_all[is.finite(r_all)]
  if (length(r_all) == 0) return(NULL)
  
  if (!is.finite(K) || K < 2) return(NULL)
  
  method <- tolower(method)
  
  if (method == "quantile") {
    probs <- seq(0, 1, length.out = K + 1)
    brks <- as.numeric(stats::quantile(r_all, probs = probs, na.rm = TRUE, type = 7))
    brks <- unique(brks)
    if (length(brks) < 3) return(NULL)
    return(brks)
  }
  
  # default: equal_width
  r_min <- min(r_all)
  r_max <- max(r_all)
  if (!is.finite(r_min) || !is.finite(r_max) || r_min == r_max) return(NULL)
  
  seq(r_min, r_max, length.out = K + 1)
}

calc_ws_spatial <- function(radius_ref, x_ref, radius_target, x_target,
                            K = 30, method = "equal_width", min_n_per_bin = 10) {
  
  rr <- safe_as_numeric(radius_ref)
  xr <- safe_as_numeric(x_ref)
  rt <- safe_as_numeric(radius_target)
  xt <- safe_as_numeric(x_target)
  
  ok_r <- is.finite(rr) & is.finite(xr)
  ok_t <- is.finite(rt) & is.finite(xt)
  
  rr <- rr[ok_r]; xr <- xr[ok_r]
  rt <- rt[ok_t]; xt <- xt[ok_t]
  
  if (length(rr) == 0 || length(rt) == 0) return(NA_real_)
  
  brks <- make_radius_bins(c(rr, rt), K = K, method = method)
  if (is.null(brks)) return(NA_real_)
  
  # 동일 breaks로 bin id 부여
  bin_r <- cut(rr, breaks = brks, include.lowest = TRUE, right = TRUE, labels = FALSE)
  bin_t <- cut(rt, breaks = brks, include.lowest = TRUE, right = TRUE, labels = FALSE)
  
  # bin별 mean + n
  dt_r <- data.table::data.table(bin = bin_r, x = xr)[is.finite(bin)]
  dt_t <- data.table::data.table(bin = bin_t, x = xt)[is.finite(bin)]
  
  prof_r <- dt_r[, .(mean_ref_bin = safe_mean(x), n_ref_bin = .N), by = bin]
  prof_t <- dt_t[, .(mean_target_bin = safe_mean(x), n_target_bin = .N), by = bin]
  
  prof <- merge(prof_r, prof_t, by = "bin", all = TRUE)
  
  # 결측/표본 부족 bin 제거 (둘 중 하나라도 부족하면 shape 비교 의미 없음)
  min_n <- as.integer(min_n_per_bin)
  if (!is.finite(min_n) || min_n < 1) min_n <- 1L
  
  prof <- prof[
    is.finite(mean_ref_bin) & is.finite(mean_target_bin) &
      (n_ref_bin >= min_n) & (n_target_bin >= min_n)
  ]
  
  if (nrow(prof) == 0) return(NA_real_)
  
  # L1 distance averaged over bins
  mean(abs(prof$mean_target_bin - prof$mean_ref_bin), na.rm = TRUE)
}

# -------------------------------------------------------------------------
# 5) Build result row (results.csv rev1 Final)
# -------------------------------------------------------------------------

make_result_row <- function(msr_name,
                            direction,
                            sigma_score,
                            cliffs_delta,
                            ws_spatial,
                            mean_ref, sd_ref, mean_target, sd_target) {
  
  data.table::data.table(
    msr         = as.character(msr_name),
    direction   = as.character(direction),
    sigma_score = as.numeric(sigma_score),
    cliffs_delta= as.numeric(cliffs_delta),
    ws_spatial  = as.numeric(ws_spatial),
    mean_ref    = as.numeric(mean_ref),
    sd_ref      = as.numeric(sd_ref),
    mean_target = as.numeric(mean_target),
    sd_target   = as.numeric(sd_target)
  )
}
