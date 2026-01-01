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
# - spatial_drift    : Wasserstein Distance (Geometry / Optimal Transport)
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
#    - Mann–Whitney U 기반 Cliff's delta (Target > Ref 우세면 +)
#    - delta = 2U/(n_t*n_r) - 1
#
#    주의:
#    - R의 wilcox.test(x, y)$statistic 은 이름이 "W"로 나오지만,
#      2-sample에서는 R 내부 정의상 U(= rank-sum에서 m(m+1)/2 뺀 값)에 해당하는 형태로 반환됨.
#    - 따라서 여기서는 추가로 m(m+1)/2를 빼지 않는다.
#    - Target 우세를 +로 만들려면 wilcox.test(x_target, x_ref)로 호출.
# -------------------------------------------------------------------------

calc_cliffs_delta <- function(x_ref, x_target) {
  xr <- suppressWarnings(as.numeric(x_ref))
  xt <- suppressWarnings(as.numeric(x_target))
  
  xr <- xr[is.finite(xr)]
  xt <- xt[is.finite(xt)]
  
  n_r <- length(xr)
  n_t <- length(xt)
  
  if (n_r < 1 || n_t < 1) return(NA_real_)
  
  # R의 "W" (2-sample): 사실상 U 형태(shifted rank-sum)로 리턴되는 값
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
# 4) spatial_drift (Sinkhorn OT on GROUP-averaged wafer maps)
#    - 목적: REF vs TARGET의 "형상(Shape)" 차이를 물리적 이동 비용으로 정량화
#
#    Phase 1) clipping (상위 20%만 남김, 나머지는 0)
#       f_clip(v) = max(v - tau, 0),  tau = quantile(v, q_clip)
#
#    Phase 2) density normalization (총합=1)
#       sum(w)==0 이면 uniform 분포 할당
#
#    Phase 3) Sinkhorn distance (entropic OT)
#       cost C_ij = ||s_i - t_j||^2  (scaled to avoid numeric underflow)
#
#    params (run.R 권장):
#      q_clip     : OT_Q           (default 0.8)
#      epsilon    : OT_EPSILON     (default 0.1)  # cost scale 적용 전제
#      max_iter   : OT_MAX_ITER    (<100)
#      cost_scale : OT_COST_SCALE  (NULL이면 자동)
#      tiny       : OT_TINY
# -------------------------------------------------------------------------

clip_relu <- function(v, q = 0.8) {
  if (length(v) == 0) return(numeric(0))
  tau <- stats::quantile(v, probs = q, na.rm = TRUE, names = FALSE, type = 7)
  pmax(v - tau, 0)
}

normalize_prob <- function(w, tiny = 1e-12) {
  n <- length(w)
  if (n == 0) return(numeric(0))
  s <- sum(w)
  if (!is.finite(s) || s <= tiny) {
    return(rep(1 / n, n))
  }
  w / s
}

sinkhorn_cost_rect <- function(p, q, C, epsilon = 0.1, max_iter = 80, tiny = 1e-12) {
  # p: length n, q: length m, C: n x m
  K <- exp(-C / epsilon)
  K[K < tiny] <- tiny
  
  u <- rep(1, length(p))
  v <- rep(1, length(q))
  
  for (t in seq_len(max_iter)) {
    Kv <- K %*% v
    Kv[ Kv < tiny ] <- tiny
    u <- p / as.numeric(Kv)
    
    KTu <- t(K) %*% u
    KTu[ KTu < tiny ] <- tiny
    v <- q / as.numeric(KTu)
  }
  
  # transport plan P_ij = u_i * K_ij * v_j
  # cost = sum(P * C)
  P <- (u * (K * rep(v, each = length(u))))
  sum(P * C)
}

# group-averaged map (칩 좌표별 평균)
make_group_mean_map <- function(dt, msr_col, group_col, group_label, x_col, y_col) {
  dt_sub <- dt[get(group_col) == group_label, .(
    x = suppressWarnings(as.numeric(get(x_col))),
    y = suppressWarnings(as.numeric(get(y_col))),
    v = suppressWarnings(as.numeric(get(msr_col)))
  )]
  
  dt_sub <- dt_sub[is.finite(x) & is.finite(y) & is.finite(v)]
  if (nrow(dt_sub) == 0) {
    return(data.table::data.table(x = numeric(0), y = numeric(0), v = numeric(0)))
  }
  
  # 좌표별 평균 (칩 pooling -> wafer 평균맵 1장)
  dt_sub[, .(v = mean(v)), by = .(x, y)]
}

calc_spatial_drift_sinkhorn <- function(map_ref, map_target,
                                        q_clip = 0.8,
                                        epsilon = 0.1,
                                        max_iter = 80,
                                        cost_scale = NULL,
                                        tiny = 1e-12) {
  # map_ref / map_target: data.table(x,y,v)
  
  if (nrow(map_ref) == 0 || nrow(map_target) == 0) return(NA_real_)
  
  # Phase 1: clipping (각 그룹 내부에서 threshold)
  w_ref <- clip_relu(map_ref$v, q = q_clip)
  w_tar <- clip_relu(map_target$v, q = q_clip)
  
  # Phase 2: normalize
  p <- normalize_prob(w_ref, tiny = tiny)
  q <- normalize_prob(w_tar, tiny = tiny)
  
  # support coords
  Xr <- map_ref$x; Yr <- map_ref$y
  Xt <- map_target$x; Yt <- map_target$y
  
  # cost scale: 좌표 스케일을 1 근처로 맞춰서 exp(-C/eps) 언더플로우 방지
  if (is.null(cost_scale)) {
    rx <- max(c(Xr, Xt)) - min(c(Xr, Xt))
    ry <- max(c(Yr, Yt)) - min(c(Yr, Yt))
    cost_scale <- max(rx, ry)
    if (!is.finite(cost_scale) || cost_scale <= tiny) cost_scale <- 1
  }
  
  dx <- outer(Xr, Xt, "-") / cost_scale
  dy <- outer(Yr, Yt, "-") / cost_scale
  C  <- dx*dx + dy*dy
  
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
