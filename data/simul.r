# data/임시_full.R
# ------------------------------------------------------------
# Fixed-R wafer, integer grid chips inside a disk.
# - You set WAFER_RADIUS first.
# - Chips are ALL integer points (X,Y) inside disk:
#     (X-CX)^2 + (Y-CY)^2 <= R^2
#   with X,Y in [0, 2R] and center (CX,CY) = (R,R).
# - All ROOTIDs share the exact same (X,Y,Chip) layout (NO sampling).
# - Outputs:
#   * data/ROOTID.csv    : ROOTID, GROUP
#   * data/raw.csv       : ROOTID, X, Y, Chip, Radius, PARTID, MSR000..MSR299
#   * data/truth_msr.csv : MSR answer key (effect type + params)
# ------------------------------------------------------------

suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", type = "binary")
  }
})
library(data.table)

set.seed(42)

# -----------------------------
# knobs (YOU set R first)
# -----------------------------
WAFER_RADIUS <- 17
CX <- WAFER_RADIUS
CY <- WAFER_RADIUS

N_LOT           <- 3
N_WAFER_PER_LOT <- 10

N_MSR <- 300
LOT_PREFIX    <- "AAA"
LOT_START_NUM <- 0

p_B <- 0.5  # wafer-level GROUP B ratio

# MSR effect mix (per MSR)
p_effect <- c(none=0.55, mean=0.15, scale=0.15, shape=0.10, hotspot=0.05)

# effect magnitudes (global defaults; logged in truth file)
delta_mean   <- 0.35
scale_mult   <- 1.35
mix_prob     <- 0.12
hot_strength <- 1.8

# -----------------------------
# helpers
# -----------------------------
points_in_disk_grid <- function(R, cx=R, cy=R) {
  xs <- 0:(2*R)
  ys <- 0:(2*R)
  grid <- CJ(X = xs, Y = ys)
  grid[(X - cx)^2 + (Y - cy)^2 <= R^2]
}

smooth_field <- function(X, Y, n_bumps = 6, sigma_bump = 6, amp_sd = 1.0,
                         cx, cy, R) {
  centers_x <- runif(n_bumps, cx - 0.7*R, cx + 0.7*R)
  centers_y <- runif(n_bumps, cy - 0.7*R, cy + 0.7*R)
  amps <- rnorm(n_bumps, mean = 0, sd = amp_sd)
  
  f <- rep(0, length(X))
  for (j in seq_len(n_bumps)) {
    d2 <- (X - centers_x[j])^2 + (Y - centers_y[j])^2
    f <- f + amps[j] * exp(-d2 / (2 * sigma_bump^2))
  }
  f
}

base_msr_spatial <- function(dt_xy, noise_sd = 0.30, heavy_tail_prob = 0.08,
                             cx, cy, R) {
  X <- dt_xy$X
  Y <- dt_xy$Y
  
  f <- smooth_field(
    X, Y,
    n_bumps    = sample(4:8, 1),
    sigma_bump = sample(c(4, 5, 6, 8, 10), 1),
    amp_sd     = runif(1, 0.8, 1.4),
    cx = cx, cy = cy, R = R
  )
  
  u <- runif(length(X))
  eps <- ifelse(
    u < (1 - heavy_tail_prob),
    rnorm(length(X), 0, noise_sd),
    rt(length(X), df=3) * noise_sd * 2.0
  )
  
  f + eps
}

apply_effect_spatial <- function(x, dt_xy, group, effect_type,
                                 cx, cy, R,
                                 delta_mean, scale_mult, mix_prob, hot_strength) {
  gB <- (group == "B")
  if (effect_type == "none") return(x)
  
  if (effect_type == "mean") {
    x[gB] <- x[gB] + delta_mean
    return(x)
  }
  
  if (effect_type == "scale") {
    x[gB] <- x[gB] * scale_mult
    return(x)
  }
  
  if (effect_type == "shape") {
    idx <- which(gB)
    u <- runif(length(idx))
    add <- ifelse(u < mix_prob,
                  rnorm(length(idx), 0, 1.5),
                  rnorm(length(idx), 0, 0.2))
    x[idx] <- x[idx] + add
    return(x)
  }
  
  if (effect_type == "hotspot") {
    X <- dt_xy$X
    Y <- dt_xy$Y
    
    shift <- max(2, round(0.35 * R))
    cA <- c(cx + shift, cy + round(0.15 * R))
    cB <- c(cx - shift, cy - round(0.15 * R))
    
    dA <- (X - cA[1])^2 + (Y - cA[2])^2
    dB <- (X - cB[1])^2 + (Y - cB[2])^2
    
    sig <- max(3, round(0.20 * R))
    bumpA <- exp(-dA / (2*(sig^2))) * hot_strength
    bumpB <- exp(-dB / (2*(sig^2))) * hot_strength
    
    x[!gB] <- x[!gB] + bumpA[!gB]
    x[gB]  <- x[gB]  + bumpB[gB]
    return(x)
  }
  
  x
}

# -----------------------------
# 0) output directory
# -----------------------------
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) build common chip layout (IDENTICAL for all wafers)
# -----------------------------
dt_pts <- points_in_disk_grid(WAFER_RADIUS, CX, CY)
dt_pts[, Radius := sqrt((X - CX)^2 + (Y - CY)^2)]
dt_pts[, Chip   := paste0(X, "_", Y)]

N_CHIP_PER_WF <- nrow(dt_pts)

cat(sprintf("[Chip Layout] R=%d | chips per wafer (integer points in disk) = %d\n",
            WAFER_RADIUS, N_CHIP_PER_WF))

# -----------------------------
# 2) ROOTID meta (wafer-level)
# -----------------------------
lot_ids <- sprintf("%s%03d", LOT_PREFIX, LOT_START_NUM + 0:(N_LOT-1))
wf_nums <- sprintf("%02d", 0:(N_WAFER_PER_LOT-1))
ROOTID  <- as.vector(outer(lot_ids, wf_nums, FUN = function(l,w) paste0(l, ".", w)))

dt_meta <- data.table(
  ROOTID = ROOTID,
  GROUP  = ifelse(runif(length(ROOTID)) < p_B, "B", "A")
)

# -----------------------------
# 3) raw skeleton = (ROOTID x common chip layout)
#    - NO sampling
# -----------------------------
dt_meta[, dummy := 1L]
dt_pts[,  dummy := 1L]
dt_raw <- merge(dt_meta, dt_pts, by = "dummy", allow.cartesian = TRUE)
dt_raw[, dummy := NULL]
dt_meta[, dummy := NULL]
dt_pts[,  dummy := NULL]

# add PARTID (chip-level)
dt_raw[, PARTID := sample(c("PART_A","PART_B","PART_C"), .N, replace = TRUE)]

# reorder base columns
setcolorder(dt_raw, c("ROOTID","GROUP","X","Y","Chip","Radius","PARTID"))

# -----------------------------
# 4) MSR truth map (answer key)
# -----------------------------
msr_names    <- sprintf("MSR%03d", 0:(N_MSR-1))
effect_types <- sample(names(p_effect), N_MSR, replace = TRUE, prob = p_effect)

dt_truth <- data.table(
  MSR = msr_names,
  effect_type = effect_types,
  truth_mean_shift   = effect_types == "mean",
  truth_scale_change = effect_types == "scale",
  truth_shape_change = effect_types == "shape",
  truth_hotspot_move = effect_types == "hotspot",
  delta_mean   = ifelse(effect_types == "mean",   delta_mean,   NA_real_),
  scale_mult   = ifelse(effect_types == "scale",  scale_mult,   NA_real_),
  mix_prob     = ifelse(effect_types == "shape",  mix_prob,     NA_real_),
  hot_strength = ifelse(effect_types == "hotspot",hot_strength, NA_real_)
)

# -----------------------------
# 5) add MSR columns
#    - base centered to mean 0 so "mean" effect cleanly controls wafer-average shift
# -----------------------------
dt_xy_all <- dt_raw[, .(X, Y)]

for (k in seq_along(msr_names)) {
  nm <- msr_names[k]
  et <- effect_types[k]
  
  x <- base_msr_spatial(
    dt_xy_all,
    noise_sd = runif(1, 0.20, 0.45),
    heavy_tail_prob = runif(1, 0.03, 0.10),
    cx = CX, cy = CY, R = WAFER_RADIUS
  )
  
  x <- x - mean(x)  # center
  
  x <- apply_effect_spatial(
    x, dt_xy_all, dt_raw$GROUP, et,
    cx = CX, cy = CY, R = WAFER_RADIUS,
    delta_mean = delta_mean,
    scale_mult = scale_mult,
    mix_prob = mix_prob,
    hot_strength = hot_strength
  )
  
  dt_raw[, (nm) := round(x, 6)]
}

# raw에는 GROUP 안 남기고 싶다 했으니 제거(메타에서 관리)
dt_raw[, GROUP := NULL]

# force final column order
dt_raw <- dt_raw[, c("ROOTID","X","Y","Chip","Radius","PARTID", msr_names), with = FALSE]

# -----------------------------
# 6) write files
# -----------------------------
fwrite(dt_meta,  file.path("data", "ROOTID.csv"))
fwrite(dt_raw,   file.path("data", "raw.csv"))
fwrite(dt_truth, file.path("data", "truth_msr.csv"))

cat("Saved:\n")
cat("- data/ROOTID.csv    (", nrow(dt_meta),  " wafer rows)\n", sep="")
cat("- data/raw.csv       (", nrow(dt_raw),   " chip rows; chips/wafer=", N_CHIP_PER_WF, "; MSRs=", N_MSR, ")\n", sep="")
cat("- data/truth_msr.csv (", nrow(dt_truth), " MSR keys)\n", sep="")
cat("\nMSR effect mix:\n")
print(table(effect_types))











# 시각화 ---------------------------------------------------------------------



library(data.table)
library(ggplot2)

# -----------------------------
# 0) load
# -----------------------------
dt_raw   <- fread("data/raw.csv")
dt_meta  <- fread("data/ROOTID.csv")
dt_truth <- fread("data/truth_msr.csv")

# raw에 GROUP 붙이기
dt <- merge(dt_raw, dt_meta, by = "ROOTID", all.x = TRUE)

# -----------------------------
# 1) 예시 MSR 자동 선택 (mean / scale / hotspot / none)
# -----------------------------
pick_one <- function(effect) {
  x <- dt_truth[effect_type == effect, MSR]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

msr_examples <- na.omit(c(
  pick_one("mean"),
  pick_one("scale"),
  pick_one("hotspot"),
  pick_one("none")
))

print(msr_examples)

# -----------------------------
# 2) plotting function: GROUP별 (X,Y) 평균맵
# -----------------------------
plot_msr_ab_meanmap <- function(msr_name, dt, R = WAFER_RADIUS) {
  # (X,Y)별 그룹 평균
  dtm <- dt[, .(val = mean(get(msr_name), na.rm = TRUE)), by = .(GROUP, X, Y)]
  
  # 원 경계선(참고용): 중심이 (R,R)인 원
  theta <- seq(0, 2*pi, length.out = 360)
  circle <- data.table(
    X = R + R * cos(theta),
    Y = R + R * sin(theta)
  )
  
  ggplot(dtm, aes(x = X, y = Y, fill = val)) +
    geom_tile() +
    geom_path(data = circle, aes(x = X, y = Y), inherit.aes = FALSE) +
    coord_fixed() +
    facet_wrap(~ GROUP, nrow = 1) +
    labs(title = paste0(msr_name, " : GROUP A/B mean wafer map"),
         x = "X", y = "Y", fill = "mean") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal(base_size = 12)
}

# -----------------------------
# 3) draw (print) + optional save
# -----------------------------
for (m in msr_examples) {
  p <- plot_msr_ab_meanmap(m, dt, R = WAFER_RADIUS)
  print(p)
  
  # 파일로 저장하고 싶으면 주석 해제
  # ggsave(filename = paste0("data/meanmap_", m, ".png"),
  #        plot = p, width = 10, height = 4, dpi = 150)
}

