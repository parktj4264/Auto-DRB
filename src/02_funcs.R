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
# - spatial_drift    : Sinkhorn OT distance on preprocessed wafer maps
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

# MAP Average Function
make_group_mean_map <- function(dt, msr_col, group_col, group_label,
                                x_col = "X", y_col = "Y") {
  
  # dt: 01_load.R 결과 (GROUP 붙어있어야 함)
  # return: data.table(x, y, v)
  
  d0 <- dt[get(group_col) == group_label, .(
    x = safe_as_numeric(get(x_col)),
    y = safe_as_numeric(get(y_col)),
    v = safe_as_numeric(get(msr_col))
  )]
  
  d0 <- d0[is.finite(x) & is.finite(y) & is.finite(v)]
  if (nrow(d0) == 0) return(data.table::data.table(x=numeric(0), y=numeric(0), v=numeric(0)))
  
  # "그룹 내 웨이퍼들"을 평균내서 대표 맵 생성
  d0[, .(v = mean(v, na.rm = TRUE)), by = .(x, y)]
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
  if (!is.finite(sigma_score) || !is.finite(sigma_level)) return("Stable")
  if (sigma_score >= sigma_level) return("Up")
  if (sigma_score <= -sigma_level) return("Down")
  "Stable"
}

# -------------------------------------------------------------------------
# 3) cliffs_delta (Stochastic Dominance)
#    - delta = 2U/(n_t*n_r) - 1  (Target > Ref 우세면 +)
# -------------------------------------------------------------------------

calc_cliffs_delta <- function(x_ref, x_target) {
  xr <- suppressWarnings(as.numeric(x_ref))
  xt <- suppressWarnings(as.numeric(x_target))
  
  xr <- xr[is.finite(xr)]
  xt <- xt[is.finite(xt)]
  
  n_r <- length(xr)
  n_t <- length(xt)
  
  if (n_r < 1 || n_t < 1) return(NA_real_)
  
  U <- tryCatch({
    as.numeric(stats::wilcox.test(xt, xr, alternative = "two.sided", exact = FALSE)$statistic)
  }, error = function(e) NA_real_)
  
  if (!is.finite(U)) return(NA_real_)
  
  denom <- n_t * n_r
  if (denom <= 0) return(NA_real_)
  
  delta <- (2 * U / denom) - 1
  as.numeric(delta)
}

# -------------------------------------------------------------------------
# 4) spatial_drift (Preprocess -> Sinkhorn OT)
#    - map_ref / map_target: data.table(x, y, v)
#    - Output: scalar distance (shape drift)
#
# Preprocess (rev1.1):
#  (1) Z-filter: abs(z) > sigma_thresh -> energy w0
#  (2) Dual convolution:
#      - sm = Smooth(w0)          : intensity diffusion
#      - d  = Smooth(mask(w0>0))  : density estimation
#  (3) final = sm * Sat(d * mask_gain)  (cluster 강조 / 고립 노이즈 억제)
#  (4) OT는 final을 질량(w)로 사용
# -------------------------------------------------------------------------

# ---------------------------------------
# [Helper 1] Gaussian Kernel (2D)
# ---------------------------------------
make_gaussian_kernel <- function(sigma = 1.0) {
  if (!is.finite(sigma) || sigma <= 0) return(matrix(1, 1, 1))
  
  k_size <- ceiling(3 * sigma) * 2 + 1
  center <- (k_size + 1) / 2
  
  ii <- seq_len(k_size)
  d2 <- (ii - center)^2
  
  kernel <- exp(-(outer(d2, d2, "+")) / (2 * sigma^2))
  kernel / sum(kernel)
}

# ---------------------------------------
# [Helper 2] 2D Convolution (same padding)
# ---------------------------------------
conv2_same <- function(mat, kernel) {
  nr <- nrow(mat); nc <- ncol(mat)
  kr <- nrow(kernel); kc <- ncol(kernel)
  khr <- floor(kr / 2); khc <- floor(kc / 2)
  
  padded <- matrix(0, nr + 2 * khr, nc + 2 * khc)
  padded[(khr + 1):(nr + khr), (khc + 1):(nc + khc)] <- mat
  
  out <- matrix(0, nr, nc)
  
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      sub_m <- padded[i:(i + 2 * khr), j:(j + 2 * khc)]
      out[i, j] <- sum(sub_m * kernel)
    }
  }
  out
}

# ---------------------------------------
# [Helper 3] Preprocess: Z-filter -> Dual Smooth -> Saturating Mask
#   - return: data.table(x, y, w)  (w >= 0)
# ---------------------------------------
preprocess_map_robust <- function(map_dt,
                                  sigma_thresh = 3.0,
                                  smooth_sigma = 1.0,
                                  mask_gain = 2.0,
                                  max_grid = 400000,
                                  tiny = 1e-12) {
  
  if (is.null(map_dt) || nrow(map_dt) == 0) {
    return(data.table::data.table(x = numeric(0), y = numeric(0), w = numeric(0)))
  }
  
  dt0 <- data.table::as.data.table(map_dt)[, .(x, y, v)]
  dt0[, x := safe_as_numeric(x)]
  dt0[, y := safe_as_numeric(y)]
  dt0[, v := safe_as_numeric(v)]
  
  # 좌표/값 유효성
  dt0 <- dt0[is.finite(x) & is.finite(y) & is.finite(v)]
  if (nrow(dt0) == 0) {
    return(data.table::data.table(x = numeric(0), y = numeric(0), w = numeric(0)))
  }
  
  # 중복좌표 평균
  dt0 <- dt0[, .(v = mean(v, na.rm = TRUE)), by = .(x, y)]
  
  vals <- dt0$v
  mu   <- mean(vals, na.rm = TRUE)
  sdv  <- stats::sd(vals, na.rm = TRUE)
  
  # sd가 거의 0이면 신호 없음 처리
  if (!is.finite(sdv) || sdv < tiny) {
    dt0[, w := 0.0]
    return(dt0[, .(x, y, w)])
  }
  
  # Z-score
  z <- (vals - mu) / sdv
  z[!is.finite(z)] <- 0
  
  if (!is.finite(sigma_thresh) || sigma_thresh < 0) sigma_thresh <- 0
  if (!is.finite(mask_gain) || mask_gain < 0) mask_gain <- 0
  
  # 에너지 필터
  w0 <- ifelse(abs(z) > sigma_thresh, abs(z), 0)
  
  # 스무딩 비활성 옵션
  if (!is.finite(smooth_sigma) || smooth_sigma <= 0) {
    dt0[, w := as.numeric(w0)]
    dt0[w < 0, w := 0.0]
    return(dt0[, .(x, y, w)])
  }
  
  # grid 구성
  xs <- sort(unique(dt0$x))
  ys <- sort(unique(dt0$y))
  
  if (length(xs) * length(ys) > max_grid) {
    # grid 폭발 방지: smoothing 스킵
    dt0[, w := as.numeric(w0)]
    dt0[w < 0, w := 0.0]
    return(dt0[, .(x, y, w)])
  }
  
  x_to_c <- match(dt0$x, xs)
  y_to_r <- match(dt0$y, ys)
  
  mat_w0 <- matrix(0, nrow = length(ys), ncol = length(xs))
  mat_w0[cbind(y_to_r, x_to_c)] <- w0
  
  mat_mask <- (mat_w0 > 0) * 1.0
  
  kernel <- make_gaussian_kernel(smooth_sigma)
  
  # dual convolution
  sm_mat <- conv2_same(mat_w0,   kernel)  # intensity
  d_mat  <- conv2_same(mat_mask, kernel)  # density
  
  # saturating density attention
  final_mat <- sm_mat * pmin(d_mat * mask_gain, 1.0)
  
  w_final <- as.numeric(final_mat[cbind(y_to_r, x_to_c)])
  w_final[!is.finite(w_final)] <- 0
  w_final[w_final < 0] <- 0
  
  dt0[, w := w_final]
  dt0[, .(x, y, w)]
}

# ---------------------------------------
# [Helper 4] Normalize to probability
# ---------------------------------------
normalize_prob <- function(w, tiny = 1e-12) {
  if (length(w) == 0) return(numeric(0))
  s <- sum(w, na.rm = TRUE)
  if (!is.finite(s) || s < tiny) return(rep(1 / length(w), length(w)))
  w / s
}

# ---------------------------------------
# [Helper 5] Sinkhorn cost for rectangular supports
#   p: length n, q: length m, C: n x m
#   return: scalar cost
# ---------------------------------------
sinkhorn_cost_rect <- function(p, q, C,
                               epsilon = 0.1,
                               max_iter = 80,
                               tiny = 1e-12) {
  n <- length(p); m <- length(q)
  if (n == 0 || m == 0) return(NA_real_)
  if (!is.finite(epsilon) || epsilon <= 0) epsilon <- 0.1
  
  K <- exp(-C / epsilon)
  K[!is.finite(K)] <- 0
  K[K < tiny] <- 0
  
  u <- rep(1, n)
  v <- rep(1, m)
  
  for (it in seq_len(max_iter)) {
    Kv <- as.numeric(K %*% v)
    Kv[Kv < tiny] <- tiny
    u <- p / Kv
    
    Ktu <- as.numeric(t(K) %*% u)
    Ktu[Ktu < tiny] <- tiny
    v <- q / Ktu
  }
  
  KCv  <- as.numeric((K * C) %*% v)
  cost <- sum(u * KCv)
  as.numeric(cost)
}

# ---------------------------------------
# [Main] spatial drift (preprocess + sinkhorn)
# ---------------------------------------
calc_spatial_drift_sinkhorn <- function(map_ref, map_target,
                                        sigma_thresh = 3.0,
                                        smooth_sigma = 1.0,
                                        mask_gain = 2.0,
                                        epsilon = 0.1,
                                        max_iter = 80,
                                        cost_scale = NULL,
                                        empty_penalty = 1.0,
                                        tiny = 1e-12) {
  if (is.null(map_ref) || is.null(map_target)) return(NA_real_)
  if (nrow(map_ref) == 0 || nrow(map_target) == 0) return(NA_real_)
  
  # Phase 1: preprocess (dual conv + saturating mask)
  dt_r <- preprocess_map_robust(map_ref,
                                sigma_thresh = sigma_thresh,
                                smooth_sigma = smooth_sigma,
                                mask_gain    = mask_gain,
                                tiny         = tiny)
  
  dt_t <- preprocess_map_robust(map_target,
                                sigma_thresh = sigma_thresh,
                                smooth_sigma = smooth_sigma,
                                mask_gain    = mask_gain,
                                tiny         = tiny)
  
  sr <- sum(dt_r$w, na.rm = TRUE)
  st <- sum(dt_t$w, na.rm = TRUE)
  
  if (sr < tiny && st < tiny) return(0.0)
  if (sr < tiny || st < tiny) return(as.numeric(empty_penalty))
  
  dt_r2 <- dt_r[w > 0]
  dt_t2 <- dt_t[w > 0]
  
  if (nrow(dt_r2) == 0 && nrow(dt_t2) == 0) return(0.0)
  if (nrow(dt_r2) == 0 || nrow(dt_t2) == 0) return(as.numeric(empty_penalty))
  
  p <- normalize_prob(dt_r2$w, tiny = tiny)
  q <- normalize_prob(dt_t2$w, tiny = tiny)
  
  Xr <- dt_r2$x; Yr <- dt_r2$y
  Xt <- dt_t2$x; Yt <- dt_t2$y
  
  # cost scaling
  if (is.null(cost_scale)) {
    rx <- max(c(Xr, Xt), na.rm = TRUE) - min(c(Xr, Xt), na.rm = TRUE)
    ry <- max(c(Yr, Yt), na.rm = TRUE) - min(c(Yr, Yt), na.rm = TRUE)
    cost_scale <- max(rx, ry)
    if (!is.finite(cost_scale) || cost_scale < tiny) cost_scale <- 1
  } else {
    if (!is.finite(cost_scale) || cost_scale < tiny) cost_scale <- 1
  }
  
  dx <- outer(Xr, Xt, "-") / cost_scale
  dy <- outer(Yr, Yt, "-") / cost_scale
  C  <- dx * dx + dy * dy
  
  sinkhorn_cost_rect(p, q, C, epsilon = epsilon, max_iter = max_iter, tiny = tiny)
}

# -------------------------------------------------------------------------
# 5) Build result row (results.csv rev1 Final)
# -------------------------------------------------------------------------

make_result_row <- function(msr_name,
                            direction,
                            sigma_score,
                            cliffs_delta,
                            spatial_drift,
                            mean_ref, sd_ref, mean_target, sd_target) {
  
  data.table::data.table(
    msr           = as.character(msr_name),
    direction     = as.character(direction),
    sigma_score   = as.numeric(sigma_score),
    cliffs_delta  = as.numeric(cliffs_delta),
    spatial_drift = as.numeric(spatial_drift),
    mean_ref      = as.numeric(mean_ref),
    sd_ref        = as.numeric(sd_ref),
    mean_target   = as.numeric(mean_target),
    sd_target     = as.numeric(sd_target)
  )
}
