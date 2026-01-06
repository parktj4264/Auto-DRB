library(data.table)
library(ggplot2)

# ==============================================================================
# [Data Loading & Preprocessing]
# - 기본 데이터(raw), 메타데이터(ROOTID), 정답지(truth)를 로드하고 병합합니다.
# - data.table의 fread로 고속 I/O를 수행합니다.
# ==============================================================================
dt_raw   <- fread("data/raw.csv")
dt_meta  <- fread("data/ROOTID.csv")
dt_truth <- fread("data/truth_msr.csv")
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# ==============================================================================
# [Helper Function: pick_one]
# - 특정 불량 유형(effect_type)에 해당하는 MSR(계측 데이터 컬럼)을 랜덤 샘플링합니다.
# - 분석 대상 컬럼을 동적으로 선택하기 위한 유틸리티 함수입니다.
# ==============================================================================
pick_one <- function(effect) {
  x <- dt_truth[effect_type == effect, MSR]
  if (length(x) == 0) return(NA_character_)
  sample(x, size = 1)
}


msr <- pick_one("hotspot")
if (is.na(msr)) msr <- pick_one("none")
stopifnot(!is.na(msr))
cat("Using MSR =", msr, "\n")

# msr <- "MSR089"

# ==============================================================================
# [Function: make_gaussian_kernel]
# - 2D Gaussian Kernel(Filter)을 생성합니다.
# - Mathematical Logic:
#     K(x, y) = (1 / (2 * pi * sigma^2)) * exp( - (x^2 + y^2) / (2 * sigma^2) )
#     (코드에서는 정규화 상수 대신 합이 1이 되도록 마지막에 나눕니다.)
# - Parameters:
#     sigma: 표준편차 (커널의 퍼짐 정도 결정)
#     k_size: 커널의 크기 (일반적으로 3*sigma 규칙을 따라 정보 손실 최소화)
# ==============================================================================
make_gaussian_kernel <- function(sigma = 1.0) {
  if (!is.finite(sigma) || sigma <= 0) return(matrix(1, 1, 1))
  
  # Kernel Size 설정: 3-sigma rule 적용 (99.7% 정보 포함 범위)
  # Size = ceiling(3 * sigma) * 2 + 1 (홀수 크기 보장)
  k_size <- ceiling(3 * sigma) * 2 + 1
  center <- (k_size + 1) / 2
  ii <- seq_len(k_size)
  
  # 거리 제곱 계산: d^2 = (x - center)^2
  d2 <- (ii - center)^2
  
  # Outer Product로 2D 거리 행렬 생성 및 Gaussian 함수 적용
  # kernel ~ exp( - (dx^2 + dy^2) / (2 * sigma^2) )
  kernel <- exp(-(outer(d2, d2, "+")) / (2 * sigma^2))
  
  # Normalize: 커널의 합이 1이 되도록 조정 (Brightness 보존)
  # sum(K) = 1
  kernel / sum(kernel)
}

# ==============================================================================
# [Function: conv2_same]
# - 2D Convolution (합성곱) 연산을 수행합니다.
# - Output size가 Input size와 동일하도록 'Same Padding'을 적용합니다.
# - Mathematical Logic:
#     (f * g)[x, y] = sum_{i, j} ( f[x-i, y-j] * g[i, j] )
# - Implementation:
#     1. Zero-padding 수행 (가장자리 처리)
#     2. Sliding Window 방식으로 커널과 부분 행렬의 내적(Dot Product) 합 계산
# ==============================================================================
conv2_same <- function(mat, kernel) {
  nr <- nrow(mat); nc <- ncol(mat)
  kr <- nrow(kernel); kc <- ncol(kernel)
  khr <- floor(kr / 2); khc <- floor(kc / 2)
  
  # Zero-padding: 커널의 반경만큼 상하좌우에 0을 채움
  padded <- matrix(0, nr + 2*khr, nc + 2*khc)
  padded[(khr+1):(nr+khr), (khc+1):(nc+khc)] <- mat
  
  out <- matrix(0, nr, nc)
  # Convolution Loop
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      sub_m <- padded[i:(i+2*khr), j:(j+2*khc)]
      out[i, j] <- sum(sub_m * kernel)
    }
  }
  out
}

# ==============================================================================
# [Function: calc_smoothing_stages_with_raw]
# - Raw 데이터로부터 노이즈를 제거하고 주요 신호(Hotspot)를 스무딩하는 핵심 파이프라인.
# - Logic Flow:
#     1. Grid Aggregation: (x, y) 좌표별 평균값 계산
#     2. Z-Score Normalization: Z = (X - mu) / sigma
#     3. Thresholding (w0): |Z| > threshold 인 값만 남김 (High Energy Signal 추출)
#     4. Mask Generation: 신호가 있는 곳을 1, 없는 곳을 0으로 마스킹
#     5. Dual Convolution:
#        - sm_mat = w0 * Kernel (신호 자체의 스무딩, Intensity Diffusion)
#        - d_mat  = Mask * Kernel (신호 밀도/분포의 스무딩, Density Estimation)
#     6. Final Computation:
#        - Final = Smooth(Signal) * Saturation(Smooth(Mask) * Gain)
#        - 주변에 신호가 밀집된 곳(d_mat 높음)의 가중치를 살리는 Attention 메커니즘 유사
# ==============================================================================
calc_smoothing_stages_with_raw <- function(map_dt,
                                           sigma_thresh = 3.0,
                                           smooth_sigma = 1.5,
                                           # iso_alpha = 2.0,
                                           mask_gain = 3.0,
                                           tiny = 1e-12) {
  # 1. Prepare Data Grid
  dt0 <- as.data.table(map_dt)[, .(x, y, v)]
  dt0 <- dt0[is.finite(x) & is.finite(y)]
  dt0 <- dt0[, .(v = mean(v, na.rm=TRUE)), by=.(x,y)]
  
  # 2. Global Statistics (Mean, SD)
  vals <- dt0$v
  mu <- mean(vals, na.rm=TRUE)
  sdv <- sd(vals, na.rm=TRUE)
  
  # 예외 처리: 표준편차가 거의 0인 경우 (모든 값이 동일)
  if (!is.finite(sdv) || sdv < tiny) {
    dt0[, `:=`(raw=v, w0=0, sm=0, d=0, final=0)]
  } else {
    # 3. Z-score Calculation
    # Z = (v - mu) / sdv
    z <- (vals - mu) / sdv
    z[!is.finite(z)] <- 0
    
    if (!is.finite(sigma_thresh) || sigma_thresh < 0) sigma_thresh <- 0
    
    # 4. Energy Filtering (w0)
    # 배경 노이즈 제거: 특정 Sigma 이상 튀는 값(Outlier)만 남김
    w0 <- ifelse(abs(z) > sigma_thresh, abs(z), 0)
    
    # 5. Coordinate Mapping (Sparse to Dense Matrix)
    xs <- sort(unique(dt0$x))
    ys <- sort(unique(dt0$y))
    x_to_c <- match(dt0$x, xs)
    y_to_r <- match(dt0$y, ys)
    
    mat_w0 <- matrix(0, nrow=length(ys), ncol=length(xs))
    mat_w0[cbind(y_to_r, x_to_c)] <- w0
    
    # Mask Matrix: 신호 유무 (Binary)
    mat_mask <- (mat_w0 > 0) * 1.0
    
    # 6. Spatial Smoothing (Convolution)
    kernel <- make_gaussian_kernel(smooth_sigma)
    
    # sm_mat: 신호 강도의 확산 (Signal Diffusion)
    sm_mat <- conv2_same(mat_w0, kernel)
    
    # d_mat: 신호 밀도의 확산 (Density Estimation)
    d_mat  <- conv2_same(mat_mask, kernel)
    
    # 7. Final Composite Score
    # final = sm_mat * min(d_mat * mask_gain, 1.0)
    # 해석: 스무딩된 신호 강도에 "주변 밀도 포화 함수"를 곱하여, 
    #       고립된 노이즈는 억제하고 군집된 불량(Cluster)을 강조함.
    # # final_mat <- sm_mat * (d_mat ^ iso_alpha)
    final_mat <- sm_mat * pmin(d_mat * mask_gain, 1.0)
    
    
    # 결과 매핑: Matrix -> Data Table
    dt0[, raw := v]
    dt0[, w0 := w0]
    dt0[, sm := as.numeric(sm_mat[cbind(y_to_r, x_to_c)])]
    dt0[, d  := as.numeric(d_mat[cbind(y_to_r, x_to_c)])]
    dt0[, final := as.numeric(final_mat[cbind(y_to_r, x_to_c)])]
  }
  
  # 8. Melting for Visualization
  out <- melt(
    dt0,
    id.vars=c("x","y"),
    measure.vars=c("raw","w0","sm","d","final"),
    variable.name="stage",
    value.name="val"
  )
  
  # 9. Post-processing & Scaling (Min-Max Normalization per stage)
  # Scaled = (val - min) / (max - min)
  out[, val := fifelse(!is.finite(val), 0, val)]
  out[, val_scaled := {
    mn <- min(val, na.rm=TRUE)
    mx <- max(val, na.rm=TRUE)
    if (!is.finite(mx - mn) || (mx - mn) < tiny) rep(0, .N) else (val - mn) / (mx - mn)
  }, by=stage]
  
  out[]
}

# ==============================================================================
# [Execution Block]
# - 선택된 MSR에 대해 GROUP 별로(예: Wafer, Lot) 평균 맵을 생성하고
# - 위에서 정의한 스무딩 파이프라인을 실행합니다.
# ==============================================================================

# GROUP mean map for chosen MSR
dt_map <- dt[, .(v = mean(get(msr), na.rm=TRUE)), by=.(GROUP, X, Y)]
setnames(dt_map, c("X","Y"), c("x","y"))

# Hyperparameters
sigma_thresh <- 1.0  # Z-score Cutoff (1-sigma 이상만 신호로 간주)
smooth_sigma <- 1.5  # Gaussian Kernel Sigma (스무딩 반경)
# iso_alpha    <- 1.5
mask_gain    <- 2    # Density Sensitivity Gain


stages <- dt_map[, calc_smoothing_stages_with_raw(.SD,
                                                  sigma_thresh=sigma_thresh,
                                                  smooth_sigma=smooth_sigma,
                                                  # iso_alpha=iso_alpha
                                                  mask_gain = mask_gain
),
by=.(GROUP)]

# Factor labeling for Plot Readability
stages[, stage := factor(stage,
                         levels=c("raw","w0","sm","d","final"),
                         labels=c("raw (GROUP mean map)",
                                  "w0 (Z-filter energy)",
                                  "Smooth(w0)",
                                  "d = Smooth(mask)",
                                  "final = Smooth(w0) * Sat(d)"))]

# ==============================================================================
# [Visualization]
# - ggplot2를 이용한 Facet Grid Heatmap 생성
# - viridis 컬러맵 사용 (시각적 왜곡 최소화)
# ==============================================================================
p <- ggplot(stages, aes(x=x, y=y, fill=val_scaled)) +
  geom_tile() +
  coord_fixed() +
  facet_grid(stage ~ GROUP) +
  scale_fill_viridis_c(option="C", direction=1) +
  labs(
    title = sprintf("Smoothing (Saturating Mask) | %s | sigma=%.1f, smooth=%.1f, gain=%.1f",
                    msr, sigma_thresh, smooth_sigma, mask_gain),
    x="X", y="Y", fill="scaled"
  ) +
  theme_minimal(base_size=12)

print(p)