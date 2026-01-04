library(data.table)
library(ggplot2)
library(patchwork)

# =========================================================================
# 1. [Engine] 3-Sigma Filter + Smooth (Prompt 코드 그대로)
# =========================================================================
make_gaussian_kernel <- function(sigma=1.0) {
  k_size <- ceiling(3 * sigma) * 2 + 1
  center <- (k_size + 1) / 2
  kernel <- matrix(0, k_size, k_size)
  for (i in 1:k_size) for (j in 1:k_size) kernel[i, j] <- exp(-((i-center)^2+(j-center)^2)/(2*sigma^2))
  kernel / sum(kernel)
}

preprocess_map_smooth <- function(dt, v_col="v", sigma_thresh=3.0, smooth_sigma=1.0) {
  # (네가 준 코드 그대로 사용)
  vals <- dt[[v_col]]
  if (length(vals) == 0 || sd(vals, na.rm=TRUE) == 0) return(numeric(length(vals)))
  
  z_scores <- as.numeric(scale(vals))
  # NA 방어 로직 살짝 추가 (실데이터엔 NA 있을 수 있음)
  z_scores[is.na(z_scores)] <- 0
  
  filtered <- ifelse(abs(z_scores) > sigma_thresh, abs(z_scores), 0) 
  
  # Matrix Mapping (좌표계가 정수가 아닐 수도 있으니 안전하게 매핑)
  x_range <- range(dt$x, na.rm=TRUE); y_range <- range(dt$y, na.rm=TRUE)
  # 1-based index로 변환 (정수 좌표 가정)
  idx_x <- round(dt$x - x_range[1] + 1)
  idx_y <- round(dt$y - y_range[1] + 1)
  
  rows <- max(idx_y); cols <- max(idx_x)
  mat <- matrix(0, rows, cols)
  mat[cbind(idx_y, idx_x)] <- filtered
  
  kernel <- make_gaussian_kernel(smooth_sigma)
  k_h <- floor(nrow(kernel)/2)
  padded <- matrix(0, rows+2*k_h, cols+2*k_h)
  padded[(k_h+1):(rows+k_h), (k_h+1):(cols+k_h)] <- mat
  
  res <- matrix(0, rows, cols)
  for(i in 1:rows) for(j in 1:cols) res[i,j] <- sum(padded[i:(i+2*k_h), j:(j+2*k_h)] * kernel)
  
  return(as.numeric(res[cbind(idx_y, idx_x)]))
}

# (Wrapper) 기존 함수 이름 매핑
calc_drift_pipeline <- function(...) {
  # 이전에 정의한 calc_spatial_drift_sinkhorn 함수가 있다고 가정
  calc_spatial_drift_sinkhorn(...)
}

# =========================================================================
# 2. [Data Switch] 시뮬레이션 -> 실데이터(dt) 로 교체!
# =========================================================================
# ★ 여기만 수정하면 됨!
target_msr <- "MSR190"   # 분석할 컬럼 이름
ref_group  <- "A"        # Reference 그룹명
tgt_group  <- "B"        # Target 그룹명 (불량 의심)

# (1) 데이터 추출 및 그룹별 평균 (Chip -> Wafer Map)
#     X, Y 대소문자 주의 (여기선 dt에 X, Y가 있다고 가정하고 x, y로 소문자 변환)
d_r <- dt[GROUP == ref_group, .(v = mean(get(target_msr), na.rm=TRUE)), by=.(x=X, y=Y)]
d_t <- dt[GROUP == tgt_group, .(v = mean(get(target_msr), na.rm=TRUE)), by=.(x=X, y=Y)]

# (2) 시나리오 리스트 구성 (실데이터는 1개 케이스)
scenarios <- list(
  list(name = paste("Real Data:", target_msr), ref_dt = d_r, tgt_dt = d_t)
)

# =========================================================================
# 3. 실행 및 시각화 (코드 구조 동일 유지)
# =========================================================================
plot_list <- list()

for (scn in scenarios) {
  # 데이터 준비
  d_ref <- scn$ref_dt
  d_tgt <- scn$tgt_dt
  
  # ★ 계산 (3-Sigma Filter -> Smooth)
  #    (만약 데이터가 비었으면 점수 NA 처리)
  if (nrow(d_ref) > 0 && nrow(d_tgt) > 0) {
    score <- calc_drift_pipeline(d_ref, d_tgt, sigma_thresh=3.5, smooth_sigma=13.0)
  } else {
    score <- NA
  }
  
  # 시각화 데이터 생성 (Target 기준)
  # AI View: 전처리 과정을 거친 후의 모습
  ai_view_t <- preprocess_map_smooth(d_tgt, sigma_thresh=3.5, smooth_sigma=13.0)
  
  # 실데이터의 통계량 (시각화 범위 설정용)
  mu_val <- mean(d_tgt$v, na.rm=TRUE)
  
  # [Raw Data Plot]
  # ★ coord_fixed()와 width/height=1 추가 -> 격자 깨짐 방지
  p_raw <- ggplot(d_tgt, aes(x, y, fill=v)) +
    geom_tile(width=5, height=5) + 
    coord_fixed() + 
    scale_fill_gradient2(low="blue", mid="gray90", high="red", midpoint=mu_val) +
    labs(title=sprintf("[%s] Raw Data", scn$name), 
         subtitle=sprintf("Actual Value (Mean: %.2f)", mu_val)) +
    theme_void() + theme(legend.position="bottom")
  
  # [AI View Plot]
  dt_ai <- copy(d_tgt)[, v := ai_view_t]
  p_ai <- ggplot(dt_ai, aes(x, y, fill=v)) +
    geom_tile(width=5, height=5) + 
    coord_fixed() +
    scale_fill_gradient(low="white", high="darkred") +
    labs(title=sprintf("OT Score: %.4f", score), 
         subtitle="AI View (Magnitude > 3-Sigma)") +
    theme_void() + theme(legend.position="bottom")
  
  plot_list[[1]] <- p_raw + p_ai
}

# 최종 출력
plot_list[[1]] + 
  plot_annotation(title = paste("Spatial Drift Analysis -", target_msr),
                  subtitle = "Left: Raw Group Mean / Right: What Algorithm Sees (Filtered Energy)")