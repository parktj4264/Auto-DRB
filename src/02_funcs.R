# 02_funcs.R

# src/02_funcs.R ------------------------------------------------------------

# 목적:
# - 03_calc.R에서 MSR 루프를 돌릴 때 사용할 함수 모음
# - REF/TARGET 비교 기준은 run.R의 GROUP_REF_LABEL / GROUP_TARGET_LABEL
# - 결과는 results.csv v1.0 (snake_case, ref/target 표기) 기준

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
# 2) Sigma_shift flags (k-sigma rule)
# -------------------------------------------------------------------------

sigma_shift_flags <- function(mean_ref, sd_ref, mean_target, sd_target, sigma_level) {
  # 결측/비정상 상황은 FALSE로 처리 (실무: 튼튼하게)
  if (!is.finite(mean_ref) || !is.finite(mean_target) || !is.finite(sigma_level)) {
    return(list(sigma_up = FALSE, sigma_down = FALSE))
  }
  if (!is.finite(sd_ref) || !is.finite(sd_target)) {
    return(list(sigma_up = FALSE, sigma_down = FALSE))
  }
  
  k <- sigma_level
  
  sigma_up   <- (mean_ref + k * sd_ref) < (mean_target - k * sd_target)
  sigma_down <- (mean_ref - k * sd_ref) > (mean_target + k * sd_target)
  
  list(
    sigma_up   = isTRUE(sigma_up),
    sigma_down = isTRUE(sigma_down)
  )
}

# -------------------------------------------------------------------------
# 3) Wilcoxon flag (chip pooling)
# -------------------------------------------------------------------------

wilcox_flag <- function(x_ref, x_target, alpha) {
  # NA 제거 후 표본수 너무 작으면 판단 불가 -> FALSE
  x_ref2    <- x_ref[is.finite(x_ref)]
  x_target2 <- x_target[is.finite(x_target)]
  
  if (length(x_ref2) < 2 || length(x_target2) < 2) {
    return(list(flag = FALSE, p_value = NA_real_))
  }
  
  p <- tryCatch({
    stats::wilcox.test(x_ref2, x_target2, alternative = "two.sided", exact = FALSE)$p.value
  }, error = function(e) NA_real_)
  
  list(
    flag    = is.finite(p) && (p < alpha),
    p_value = p
  )
}

# -------------------------------------------------------------------------
# 4) KS flag (Radius-weighted KS)
#    - 공간(좌/우 + 거리) 정보를 값에 섞어서 1D 분포로 만든 뒤 KS 수행
#    - z = x * (radius / max_abs_radius)
# -------------------------------------------------------------------------

make_radius_weighted_values <- function(radius, values, eps = 1e-12) {
  r <- safe_as_numeric(radius)
  x <- safe_as_numeric(values)
  
  ok <- is.finite(r) & is.finite(x)
  r <- r[ok]
  x <- x[ok]
  
  if (length(r) == 0) return(numeric(0))
  
  max_abs <- max(abs(r))
  if (!is.finite(max_abs) || max_abs < eps) return(numeric(0))
  
  # 위치 정보를 연속 가중치로 섞기: 좌(-) / 우(+) + 거리
  z <- x * (r / max_abs)
  z[is.finite(z)]
}

ks_flag_radius_weighted <- function(radius_ref, x_ref, radius_target, x_target,
                                    alpha, min_n = 30) {
  
  z_ref    <- make_radius_weighted_values(radius_ref, x_ref)
  z_target <- make_radius_weighted_values(radius_target, x_target)
  
  # 표본이 너무 작으면 안정적으로 FALSE 처리
  if (length(z_ref) < min_n || length(z_target) < min_n) {
    return(list(flag = FALSE, p_value = NA_real_, n_ref = length(z_ref), n_target = length(z_target)))
  }
  
  p <- tryCatch({
    stats::ks.test(z_ref, z_target, alternative = "two.sided")$p.value
  }, error = function(e) NA_real_)
  
  list(
    flag     = is.finite(p) && (p < alpha),
    p_value  = p,
    n_ref    = length(z_ref),
    n_target = length(z_target)
  )
}


# -------------------------------------------------------------------------
# 5) Build result row (results.csv v1.0)
# -------------------------------------------------------------------------

make_result_row <- function(msr_name,
                            sigma_up, sigma_down,
                            wilcox_flag, ks_flag,
                            mean_ref, sd_ref, mean_target, sd_target,
                            mean_diff, median_diff) {
  
  data.table::data.table(
    msr         = msr_name,
    sigma_up    = as.logical(sigma_up),
    sigma_down  = as.logical(sigma_down),
    wilcox_flag = as.logical(wilcox_flag),
    ks_flag     = as.logical(ks_flag),
    mean_ref    = as.numeric(mean_ref),
    sd_ref      = as.numeric(sd_ref),
    mean_target = as.numeric(mean_target),
    sd_target   = as.numeric(sd_target),
    mean_diff   = as.numeric(mean_diff),
    median_diff = as.numeric(median_diff)
  )
}
