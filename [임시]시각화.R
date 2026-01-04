library(data.table)
library(ggplot2)
library(patchwork)

# =========================================================================
# Part 1. [Core Engine] Math & Signal Processing Helpers
# =========================================================================

# 안전한 숫자 변환 (데이터 타입 에러 방지)
safe_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

# 1. Gaussian Kernel Generator (붓 만들기)
make_gaussian_kernel <- function(sigma = 10.0) {
  if (!is.finite(sigma) || sigma <= 0) return(matrix(1, 1, 1))
  
  k_size <- ceiling(3 * sigma) * 2 + 1  # 3-sigma rule
  center <- (k_size + 1) / 2
  
  ii <- seq_len(k_size); jj <- seq_len(k_size)
  d2 <- (ii - center)^2
  kernel <- exp(-(outer(d2, d2, "+")) / (2 * sigma^2))
  kernel / sum(kernel)
}

# 2. 2D Convolution (Naive Implementation for R)
conv2_same <- function(mat, kernel) {
  nr <- nrow(mat); nc <- ncol(mat)
  kr <- nrow(kernel); kc <- ncol(kernel)
  khr <- floor(kr / 2); khc <- floor(kc / 2)
  
  padded <- matrix(0, nr + 2 * khr, nc + 2 * khc)
  padded[(khr + 1):(nr + khr), (khc + 1):(nc + khc)] <- mat
  
  out <- matrix(0, nr, nc)
  # (Note: For very large grids > 500x500, consider using FFT or C++ extension)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      sub_m <- padded[i:(i + 2 * khr), j:(j + 2 * khc)]
      out[i, j] <- sum(sub_m * kernel)
    }
  }
  out
}

# 3. Robust Preprocessing: Z-Filter -> Smooth
#    - Noise는 죽이고(0), 불량 덩어리(Blob)만 남김
preprocess_map_robust <- function(map_dt, 
                                  sigma_thresh = 1.0, 
                                  smooth_sigma = 10.0, 
                                  tiny = 1e-12) {
  # map_dt expects columns: x, y, v
  if (is.null(map_dt) || nrow(map_dt) == 0) return(numeric(0))
  
  # 데이터 복사 및 타입 안전 변환
  vals <- safe_as_numeric(map_dt$v)
  
  # 예외: 변화가 없거나(sd=0) 데이터가 너무 적음 -> 전량 0 (Noise)
  sdv <- sd(vals, na.rm = TRUE)
  if (is.na(sdv) || sdv < tiny) return(numeric(length(vals)))
  
  # Step 1: Z-Score Filtering (Magnitude Extraction)
  z <- (vals - mean(vals, na.rm=TRUE)) / sdv
  z[is.na(z)] <- 0
  
  # 3-Sigma 밖의 값만 에너지로 인정 (방향 무시, 절대값)
  w0 <- ifelse(abs(z) > sigma_thresh, abs(z), 0)
  
  # Step 2: Gaussian Smoothing
  if (smooth_sigma > 0) {
    # 좌표 그리드 매핑
    xs <- map_dt$x; ys <- map_dt$y
    x_min <- min(xs); y_min <- min(ys)
    
    # 좌표를 1-based index로 변환 (정수 가정)
    # 실무 데이터가 실수(float) 좌표라면 rounding이나 binning 선행 필요
    idx_c <- round(xs - x_min + 1)
    idx_r <- round(ys - y_min + 1)
    
    rows <- max(idx_r); cols <- max(idx_c)
    
    # Grid가 너무 크면 스무딩 Skip (Performance 방어)
    if (rows * cols > 400000) return(w0)
    
    mat <- matrix(0, nrow = rows, ncol = cols)
    # 중복 좌표가 있다면 합산 혹은 덮어쓰기 (여기선 덮어쓰기)
    mat[cbind(idx_r, idx_c)] <- w0
    
    # Convolution
    kernel <- make_gaussian_kernel(smooth_sigma)
    res_mat <- conv2_same(mat, kernel)
    
    # 원래 위치의 값만 추출
    w_sm <- res_mat[cbind(idx_r, idx_c)]
    w_sm[w_sm < 0] <- 0 # 음수 방지
    return(as.numeric(w_sm))
  }
  
  return(as.numeric(w0))
}

# 4. Normalize Helper
normalize_prob <- function(w, tiny = 1e-12) {
  s <- sum(w, na.rm = TRUE)
  if (s < tiny) return(rep(1 / length(w), length(w))) # Fallback to uniform
  w / s
}

# 5. Sinkhorn Calculation
sinkhorn_cost_rect <- function(p, q, C, epsilon = 0.1, max_iter = 80, tiny = 1e-12) {
  if (length(p) == 0) return(NA_real_)
  
  K <- exp(-C / epsilon)
  K[K < tiny] <- 0
  
  u <- rep(1, length(p)); v <- rep(1, length(q))
  
  for (i in seq_len(max_iter)) {
    Kv <- as.numeric(K %*% v); Kv[Kv < tiny] <- tiny
    u <- p / Kv
    Ktu <- as.numeric(t(K) %*% u); Ktu[Ktu < tiny] <- tiny
    v <- q / Ktu
  }
  
  # Cost Calculation
  sum(u * ((K * C) %*% v))
}


# =========================================================================
# Part 2. [Main Function] Spatial Drift Calculator
# =========================================================================

calc_spatial_drift_sinkhorn <- function(map_ref, map_target,
                                        sigma_thresh = 3.0,
                                        smooth_sigma = 1.0,
                                        epsilon = 0.1,
                                        max_iter = 80) {
  # map_ref, map_target must have x, y, v
  if (nrow(map_ref) == 0 || nrow(map_target) == 0) return(NA_real_)
  
  # Phase 1: Robust Preprocessing
  w_ref <- preprocess_map_robust(map_ref, sigma_thresh, smooth_sigma)
  w_tar <- preprocess_map_robust(map_target, sigma_thresh, smooth_sigma)
  
  # [Optimization] Ghost Case (둘 다 노이즈만 있음) -> Pass
  sum_r <- sum(w_ref); sum_t <- sum(w_tar)
  if (sum_r < 1e-9 && sum_t < 1e-9) return(0.0)
  
  # [Penalty] Creation/Extinction (한쪽만 신호 있음) -> Fail (Max Cost)
  if (sum_r < 1e-9 || sum_t < 1e-9) return(1.0) 
  
  # Phase 2: Normalize
  p <- w_ref / sum_r
  q <- w_tar / sum_t
  
  # Phase 3: Sinkhorn
  # Sparsity Optimization (0인 값 제외하고 계산)
  idx_r <- which(p > 1e-9); idx_t <- which(q > 1e-9)
  
  p_sub <- p[idx_r]; q_sub <- q[idx_t]
  Xr <- map_ref$x[idx_r]; Yr <- map_ref$y[idx_r]
  Xt <- map_target$x[idx_t]; Yt <- map_target$y[idx_t]
  
  # Cost Scale
  cost_scale <- max(max(Xr, Xt), max(Yr, Yt))
  if (cost_scale == 0) cost_scale <- 1
  
  dx <- outer(Xr, Xt, "-") / cost_scale
  dy <- outer(Yr, Yt, "-") / cost_scale
  C  <- dx^2 + dy^2
  
  sinkhorn_cost_rect(p_sub, q_sub, C, epsilon, max_iter)
}


# =========================================================================
# Part 3. [Batch Runner] 대량 MSR 분석기
# =========================================================================

run_batch_drift_analysis <- function(dt, msr_cols, 
                                     group_col="GROUP", ref_label="A", tgt_label="B",
                                     x_col="X", y_col="Y") {
  
  # 진행 상황바 (선택사항, 없어도 됨)
  pb <- txtProgressBar(min = 0, max = length(msr_cols), style = 3)
  
  results_list <- lapply(seq_along(msr_cols), function(i) {
    col_name <- msr_cols[i]
    setTxtProgressBar(pb, i)
    
    # 1. Group Aggregation (평균 맵 생성)
    #    data.table dynamic column selection
    d_ref <- dt[get(group_col) == ref_label, 
                .(v = mean(get(col_name), na.rm=TRUE)), 
                by = c(x_col, y_col)]
    
    d_tgt <- dt[get(group_col) == tgt_label, 
                .(v = mean(get(col_name), na.rm=TRUE)), 
                by = c(x_col, y_col)]
    
    # 2. Rename columns for Engine (X -> x, Y -> y)
    if (nrow(d_ref) > 0) setnames(d_ref, c(x_col, y_col), c("x", "y"))
    if (nrow(d_tgt) > 0) setnames(d_tgt, c(x_col, y_col), c("x", "y"))
    
    # 3. Calculate Drift
    score <- NA_real_
    if (nrow(d_ref) > 0 && nrow(d_tgt) > 0) {
      score <- calc_spatial_drift_sinkhorn(
        d_ref, d_tgt, 
        sigma_thresh = 3.0, 
        smooth_sigma = 1.0
      )
    }
    
    data.table(MSR = col_name, Score = score)
  })
  
  close(pb)
  res_dt <- rbindlist(results_list)
  setorder(res_dt, -Score, na.last=TRUE)
  return(res_dt)
}


# =========================================================================
# Part 4. [Visualizer] 시각화 도구
# =========================================================================

# 1. 랭킹 차트
plot_drift_summary <- function(res_dt, top_n = 5) {
  top_data <- head(res_dt[!is.na(Score)], top_n)
  
  ggplot(top_data, aes(x = reorder(MSR, Score), y = Score)) +
    geom_col(fill = "firebrick") +
    coord_flip() +
    labs(title = "Top Spatial Drift MSRs", 
         subtitle = "Higher score means significant shape change",
         x = "Parameter", y = "Drift Score") +
    theme_minimal() +
    theme(plot.title = element_text(face="bold"))
}

# 2. 상세 비교 (Raw vs AI View)
plot_msr_detail <- function(dt, msr_name, 
                            group_col="GROUP", ref_label="A", tgt_label="B",
                            x_col="X", y_col="Y") {
  
  # 데이터 준비
  d_ref <- dt[get(group_col) == ref_label, .(x=get(x_col), y=get(y_col), v=get(msr_name))]
  d_tgt <- dt[get(group_col) == tgt_label, .(x=get(x_col), y=get(y_col), v=get(msr_name))]
  
  # 평균 맵
  map_ref <- d_ref[, .(v = mean(v, na.rm=TRUE)), by=.(x, y)]
  map_tgt <- d_tgt[, .(v = mean(v, na.rm=TRUE)), by=.(x, y)]
  
  # AI View (Energy) 계산
  ai_ref <- preprocess_map_robust(map_ref, sigma_thresh=1.0, smooth_sigma=10.0)
  ai_tgt <- preprocess_map_robust(map_tgt, sigma_thresh=1.0, smooth_sigma=10.0)
  
  # 정렬 및 결합
  setorder(map_ref, y, x); setorder(map_tgt, y, x)
  map_ref[, v_ai := ai_ref]; map_tgt[, v_ai := ai_tgt]
  
  viz_data <- rbind(
    map_ref[, .(x, y, v, v_ai, Group="Ref (A)")],
    map_tgt[, .(x, y, v, v_ai, Group="Target (B)")]
  )
  
  # ★ 시각화 옵션 강화 (격자 만들기)
  p1 <- ggplot(viz_data, aes(x, y, fill=v)) +
    # [핵심] width/height를 1로 고정해서 빈틈 메우기
    geom_tile(width = 5, height = 5) + 
    facet_wrap(~Group) +
    # [핵심] 비율 1:1 고정 (웨이퍼 찌그러짐 방지)
    coord_fixed() + 
    scale_fill_gradient2(low="blue", mid="gray95", high="red", midpoint=median(viz_data$v, na.rm=T)) +
    labs(title=paste(msr_name, "Raw Map"), subtitle="Actual Values") +
    theme_minimal() +
    theme(panel.grid = element_blank()) # 뒤에 깔리는 모눈종이 제거 (깔끔하게)
  
  p2 <- ggplot(viz_data, aes(x, y, fill=v_ai)) +
    geom_tile(width = 5, height = 5) + # 여기도 꽉 채우기
    facet_wrap(~Group) +
    coord_fixed() +
    scale_fill_gradient(low="white", high="darkred") +
    labs(title=paste(msr_name, "AI View (Filtered)"), subtitle="What OT Algorithm Sees") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  p1 / p2
}






# --- 실행 예시 ---
# MSR 컬럼 목록 자동 추출 (MSR000 ~ MSR299)
# msr_columns <- grep("^MSR", names(dt), value = TRUE)
# batch_result <- run_batch_drift_analysis(dt, msr_columns)
# 
# print(head(batch_result, 10)) # Top 10 불량 MSR 확인
# 


# 
# # 1. 데이터 로딩 (dt에 X, Y, GROUP, MSR000... 등이 있다고 가정)
# msr_cols <- grep("^MSR", names(dt), value = TRUE)
# 
# # 2. 전체 분석 실행
# batch_res <- run_batch_drift_analysis(dt, msr_cols, x_col="X", y_col="Y")
# 
# # 3. 결과 확인
# print(head(batch_res))
# plot_drift_summary(batch_res)

# 4. 1등 불량 상세 확인
top_msr <- batch_res$MSR[1]
plot_msr_detail(dt, top_msr, x_col="X", y_col="Y")



