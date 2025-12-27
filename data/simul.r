# data/임시.R
# ------------------------------------------------------------
# Simulate ROOTID.csv (wafer meta) + raw.csv (chip-level)
# - LOT: 3 (AAA000~AAA002), wafers per lot: 10 (.00~.09)
# - Chips per wafer: 800
# - X,Y: non-negative integers, full-circle wafer via center (150,150)
# - Chip: "X_Y"
# - Radius: distance from center (0~150)
# - raw.csv columns: ROOTID, X, Y, Chip, Radius, PARTID, MSR000..MSR299
# - ROOTID.csv columns: ROOTID, GROUP
# - MSR generated as smooth spatial field + noise, plus group-specific effect
# ------------------------------------------------------------

suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", type = "binary")
})
library(data.table)

set.seed(42)

# -----------------------------
# knobs
# -----------------------------
WAFER_RADIUS     <- 150
CX               <- 150
CY               <- 150

N_LOT            <- 3
N_WAFER_PER_LOT  <- 10
N_CHIP_PER_WF    <- 800
N_MSR            <- 300          # MSR000 ~ MSR299

LOT_PREFIX       <- "AAA"
LOT_START_NUM    <- 0            # AAA000부터

p_B              <- 0.5          # GROUP B 비율

# MSR별 effect mix (GROUP B 특성)
p_effect <- c(none=0.55, mean=0.15, scale=0.15, shape=0.10, hotspot=0.05)

# effect magnitudes
delta_mean   <- 0.35
scale_mult   <- 1.35
mix_prob     <- 0.12
hot_strength <- 1.8

# -----------------------------
# helpers
# -----------------------------
make_points_in_disk_int_pos <- function(n, R, cx=150, cy=150) {
  # X,Y are integers in [0, 2R] centered at (cx,cy) = (R,R) -> default 150,150
  xs <- integer(0); ys <- integer(0)
  while (length(xs) < n) {
    m <- (n - length(xs)) * 6
    x <- sample.int(2*R + 1, m, replace=TRUE) + (cx - (R + 1))  # (cx-R) ~ (cx+R)
    y <- sample.int(2*R + 1, m, replace=TRUE) + (cy - (R + 1))  # (cy-R) ~ (cy+R)
    
    keep <- ((x - cx)^2 + (y - cy)^2) <= (R^2)
    xs <- c(xs, x[keep])
    ys <- c(ys, y[keep])
  }
  data.table(X = xs[1:n], Y = ys[1:n])
}

smooth_field <- function(X, Y, n_bumps = 6, sigma_bump = 35, amp_sd = 1.0,
                         cx = 150, cy = 150, R = 150) {
  # 랜덤 범프 중심(웨이퍼 내부 위주)
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
                             cx=150, cy=150, R=150) {
  X <- dt_xy$X
  Y <- dt_xy$Y
  
  f <- smooth_field(
    X, Y,
    n_bumps = sample(4:8, 1),
    sigma_bump = sample(c(25, 30, 35, 45, 60), 1),
    amp_sd = runif(1, 0.8, 1.4),
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
                                 cx=150, cy=150, R=150) {
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
    X <- dt_xy$X; Y <- dt_xy$Y
    
    # A hotspot (cx+40, cy+20), B hotspot (cx-40, cy-20)
    cA <- c(cx + 40, cy + 20)
    cB <- c(cx - 40, cy - 20)
    
    dA <- (X - cA[1])^2 + (Y - cA[2])^2
    dB <- (X - cB[1])^2 + (Y - cB[2])^2
    
    bumpA <- exp(-dA / (2*(30^2))) * hot_strength
    bumpB <- exp(-dB / (2*(30^2))) * hot_strength
    
    x[!gB] <- x[!gB] + bumpA[!gB]
    x[gB]  <- x[gB]  + bumpB[gB]
    return(x)
  }
  
  x
}

# -----------------------------
# 1) ROOTID meta (wafer-level): ROOTID + GROUP only
# -----------------------------
lot_ids <- sprintf("%s%03d", LOT_PREFIX, LOT_START_NUM + 0:(N_LOT-1))
wf_nums <- sprintf("%02d", 0:(N_WAFER_PER_LOT-1))
ROOTID  <- as.vector(outer(lot_ids, wf_nums, FUN = function(l,w) paste0(l, ".", w)))

dt_meta <- data.table(
  ROOTID = ROOTID,
  GROUP  = ifelse(runif(length(ROOTID)) < p_B, "B", "A")
)

# -----------------------------
# 2) raw (chip-level)
# -----------------------------
raw_list <- vector("list", length = nrow(dt_meta))

for (i in seq_len(nrow(dt_meta))) {
  pts <- make_points_in_disk_int_pos(N_CHIP_PER_WF, WAFER_RADIUS, cx=CX, cy=CY)
  # pts[, Radius := round(sqrt((X - CX)^2 + (Y - CY)^2), 0)]
  pts[, Radius := (X - CX) + runif(.N, -0.5, 0.5)]
  pts[, Chip   := paste0(X, "_", Y)]
  
  raw_list[[i]] <- data.table(
    ROOTID = dt_meta$ROOTID[i],
    X      = pts$X,
    Y      = pts$Y,
    Chip   = pts$Chip,
    Radius = pts$Radius,
    PARTID = sample(c("PART_A","PART_B","PART_C"), N_CHIP_PER_WF, replace = TRUE)
  )
}

dt_raw <- rbindlist(raw_list)

# -----------------------------
# 3) add MSR columns (PARTID 다음부터 MSR 시작)
# -----------------------------
msr_names    <- sprintf("MSR%03d", 0:(N_MSR-1))
effect_types <- sample(names(p_effect), N_MSR, replace = TRUE, prob = p_effect)

# GROUP 붙여서 효과 적용
dt_raw <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# 좌표 한번만 잡아두기
dt_xy_all <- dt_raw[, .(X, Y)]

for (k in seq_along(msr_names)) {
  nm <- msr_names[k]
  et <- effect_types[k]
  
  x <- base_msr_spatial(
    dt_xy_all,
    noise_sd = runif(1, 0.20, 0.45),
    heavy_tail_prob = runif(1, 0.03, 0.10),
    cx=CX, cy=CY, R=WAFER_RADIUS
  )
  
  x <- apply_effect_spatial(x, dt_xy_all, dt_raw$GROUP, et, cx=CX, cy=CY, R=WAFER_RADIUS)
  
  dt_raw[, (nm) := round(x, 6)]
}

# GROUP은 raw에 굳이 안 남기고 싶으면 지우기 (meta에서 관리)
dt_raw[, GROUP := NULL]

# 컬럼 순서 강제
dt_raw <- dt_raw[, c("ROOTID","X","Y","Chip","Radius","PARTID", msr_names), with=FALSE]

# -----------------------------
# 4) write
# -----------------------------
dir.create("data", showWarnings = FALSE, recursive = TRUE)
fwrite(dt_meta, file.path("data", "ROOTID.csv"))
fwrite(dt_raw,  file.path("data", "raw.csv"))

cat("Saved:\n")
cat("- data/ROOTID.csv  (", nrow(dt_meta), " wafer rows)\n", sep="")
cat("- data/raw.csv     (", nrow(dt_raw), " chip rows, ", length(msr_names), " MSRs)\n", sep="")
cat(sprintf("Lots=%d, wafers/lot=%d, chips/wf=%d\n", N_LOT, N_WAFER_PER_LOT, N_CHIP_PER_WF))
cat("\nMSR effect mix:\n")
print(table(effect_types))






















# 차이 ----------------------------------------------------------------------


msr_truth <- data.table(
  MSR = sprintf("MSR%03d", 0:(N_MSR-1)),
  effect_type = effect_types
)

msr_truth[]

fwrite(
  msr_truth,
  file.path("data", "MSR_truth_정답지.csv")
)














# 시각화 ---------------------------------------------------------------------
# -------------------------------------------------------------------------
# suppressWarnings({
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", type="binary")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", type="binary")

library(data.table)
library(ggplot2)

dt_raw  <- fread(file.path("data","raw.csv"))
dt_meta <- fread(file.path("data","ROOTID.csv"))
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# ---- params
msr_pick <- c("MSR287")
msr_pick <- intersect(msr_pick, names(dt))
if (length(msr_pick) == 0) stop("msr_pick이 dt에 없음")

CX <- 150; CY <- 150; R <- 150
BIN <- 10   # <- 여기만 바꾸면 됨: 1(가장 촘촘) / 5(추천) / 10(더 매끈)

# ---- chip score (선택 MSR 평균)
dt[, chip_score := rowMeans(.SD, na.rm=TRUE), .SDcols = msr_pick]

# ---- binning: (X,Y)를 BIN 단위 타일로 묶기
# 타일 좌표는 "타일의 중심"으로 잡아두면 원형 마스크가 깔끔함
dt[, `:=`(
  Xb = (X %/% BIN) * BIN + BIN/2,
  Yb = (Y %/% BIN) * BIN + BIN/2
)]

# ---- GROUP별 타일 평균
dt_tile <- dt[, .(mean_score = mean(chip_score, na.rm=TRUE), n=.N),
              by=.(GROUP, Xb, Yb)]

# ---- 원형 마스크: 타일 중심이 원 안에 있는 것만 남김
dt_tile <- dt_tile[((Xb - CX)^2 + (Yb - CY)^2) <= R^2]

# ---- 타일 격자 전체를 "완성"해서 빈칸도 타일로 채우기(NA는 흰색으로)
xs <- seq(BIN/2, 300 + BIN/2, by=BIN)
ys <- seq(BIN/2, 300 + BIN/2, by=BIN)
grid <- CJ(Xb = xs, Yb = ys)
grid <- grid[((Xb - CX)^2 + (Yb - CY)^2) <= R^2]
grid[, tmp := 1L]

# GROUP별로 완성격자 만들고 merge
groups <- sort(unique(dt_meta$GROUP))
full_grid <- rbindlist(lapply(groups, function(g){
  tmp <- copy(grid)
  tmp[, GROUP := g]
  tmp
}))

dt_plot <- merge(full_grid[, .(GROUP, Xb, Yb)],
                 dt_tile[, .(GROUP, Xb, Yb, mean_score)],
                 by=c("GROUP","Xb","Yb"), all.x=TRUE)

# ---- 원 테두리
theta <- seq(0, 2*pi, length.out=600)
circle <- data.table(x = CX + R*cos(theta), y = CY + R*sin(theta))

# ---- plot (geom_tile)
ggplot(dt_plot, aes(x=Xb, y=Yb, fill=mean_score)) +
  geom_tile(width=BIN, height=BIN) +
  geom_path(data=circle, aes(x=x, y=y), inherit.aes=FALSE, linewidth=0.4) +
  coord_equal(xlim=c(CX-R, CX+R), ylim=c(CY-R, CY+R)) +
  facet_wrap(~GROUP) +
  theme_minimal(base_size=12) +
  labs(
    title = "Group-wise Mean Wafer Map (Tiled)",
    subtitle = paste0("BIN=", BIN, ", chip_score = mean(", paste(msr_pick, collapse=", "), ")"),
    fill = "mean score"
  )









# funny simul -------------------------------------------------------------



suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", type="binary")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", type="binary")
})
library(data.table)
library(ggplot2)

dt_raw  <- fread(file.path("data","raw.csv"))
dt_meta <- fread(file.path("data","ROOTID.csv"))
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# MSR 몇 개 골라서 "칩 스코어" 만들기 (너 마음대로)
msr_pick <- c("MSR287")
msr_pick <- intersect(msr_pick, names(dt))
dt[, chip_score := rowMeans(.SD, na.rm=TRUE), .SDcols=msr_pick]

# wafer 중심/반지름 (시뮬 기준)
CX <- 150; CY <- 150; R <- 150

# 원 테두리
theta <- seq(0, 2*pi, length.out=700)
circle <- data.table(x = CX + R*cos(theta), y = CY + R*sin(theta))








BIN <- 10

dt[, `:=`(
  Xb = (X %/% BIN) * BIN + BIN/2,
  Yb = (Y %/% BIN) * BIN + BIN/2
)]

tile <- dt[, .(mean_score = mean(chip_score, na.rm=TRUE)),
           by=.(GROUP, Xb, Yb)]
tile <- tile[((Xb-CX)^2 + (Yb-CY)^2) <= R^2]

# grid 완성
xs <- seq(BIN/2, 300 + BIN/2, by=BIN)
ys <- seq(BIN/2, 300 + BIN/2, by=BIN)
grid <- CJ(Xb=xs, Yb=ys)
grid <- grid[((Xb-CX)^2 + (Yb-CY)^2) <= R^2]
groups <- sort(unique(dt$GROUP))
full_grid <- rbindlist(lapply(groups, \(g){ tmp<-copy(grid); tmp[,GROUP:=g]; tmp }))

tile_full <- merge(full_grid, tile, by=c("GROUP","Xb","Yb"), all.x=TRUE)

# Δmap 만들기 (B-A)
Amap <- tile_full[GROUP=="A", .(Xb,Yb,A=mean_score)]
Bmap <- tile_full[GROUP=="B", .(Xb,Yb,B=mean_score)]
Dmap <- merge(Amap, Bmap, by=c("Xb","Yb"), all=TRUE)
Dmap[, delta := B - A]

# plot: A, B, Δ
pA <- ggplot(tile_full[GROUP=="A"], aes(Xb,Yb,fill=mean_score)) +
  geom_tile(width=BIN,height=BIN) +
  geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.4) +
  coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
  theme_minimal(base_size=12) +
  labs(title="Wafer Map: GROUP A", fill="score")

pB <- ggplot(tile_full[GROUP=="B"], aes(Xb,Yb,fill=mean_score)) +
  geom_tile(width=BIN,height=BIN) +
  geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.4) +
  coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
  theme_minimal(base_size=12) +
  labs(title="Wafer Map: GROUP B", fill="score")

pD <- ggplot(Dmap, aes(Xb,Yb,fill=delta)) +
  geom_tile(width=BIN,height=BIN) +
  geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.4) +
  coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
  theme_minimal(base_size=12) +
  labs(title="Δ Map = (B - A)", fill="delta")

print(pA); print(pB); print(pD)










# 반경을 BIN 단위로 뭉쳐서 프로파일
RBIN <- 5
dt[, Rb := (Radius %/% RBIN) * RBIN]

prof <- dt[, .(mean_score = mean(chip_score, na.rm=TRUE),
               sd_score   = sd(chip_score, na.rm=TRUE),
               n=.N),
           by=.(GROUP, Rb)]

ggplot(prof, aes(x=Rb, y=mean_score, group=GROUP, color=GROUP)) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.5) +
  theme_minimal(base_size=12) +
  labs(title="Radial Profile (by GROUP)", x="Radius bin", y="Mean chip_score")








wf_sum <- dt[, .(
  wf_mean = mean(chip_score, na.rm=TRUE),
  wf_sd   = sd(chip_score, na.rm=TRUE)
), by=.(ROOTID, GROUP)]

ggplot(wf_sum, aes(x=GROUP, y=wf_mean, fill=GROUP)) +
  geom_boxplot(outlier.alpha=0.3) +
  theme_minimal(base_size=12) +
  labs(title="Per-wafer mean distribution", y="wf_mean")









# 각 GROUP별 가중 중심
cm <- dt[, .(
  x_cm = sum(X * pmax(chip_score, 0), na.rm=TRUE) / sum(pmax(chip_score, 0), na.rm=TRUE),
  y_cm = sum(Y * pmax(chip_score, 0), na.rm=TRUE) / sum(pmax(chip_score, 0), na.rm=TRUE)
), by=GROUP]

ggplot() +
  geom_path(data=circle, aes(x=x,y=y), linewidth=0.5) +
  geom_point(data=cm, aes(x=x_cm, y=y_cm, color=GROUP), size=4) +
  coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
  theme_minimal(base_size=12) +
  labs(title="Center of Mass of chip_score (by GROUP)", x="X", y="Y")


























suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", type="binary")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", type="binary")
})
library(data.table)
library(ggplot2)

# -----------------------------
# load
# -----------------------------
dt_raw  <- fread(file.path("data","raw.csv"))
dt_meta <- fread(file.path("data","ROOTID.csv"))
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# wafer geometry
CX <- 150; CY <- 150; R <- 150
BIN <- 10

# circle outline
theta <- seq(0, 2*pi, length.out=700)
circle <- data.table(x = CX + R*cos(theta), y = CY + R*sin(theta))

# MSR list
msr_names <- grep("^MSR\\d{3}$", names(dt), value = TRUE)
if (length(msr_names) == 0) stop("MSR columns not found")

# precompute tile keys
dt[, `:=`(
  Xb = (X %/% BIN) * BIN + BIN/2,
  Yb = (Y %/% BIN) * BIN + BIN/2
)]

# full grid inside wafer (for consistent delta map)
xs <- seq(BIN/2, 300 + BIN/2, by=BIN)
ys <- seq(BIN/2, 300 + BIN/2, by=BIN)
grid <- CJ(Xb=xs, Yb=ys)
grid <- grid[((Xb-CX)^2 + (Yb-CY)^2) <= R^2]

# -----------------------------
# hotspot score function
# -----------------------------
hotspot_metrics <- function(tile_delta, cx=150, cy=150) {
  # tile_delta: data.table with Xb,Yb,delta (may include NA)
  d <- tile_delta[!is.na(delta)]
  
  if (nrow(d) < 10) {
    return(list(
      hotspot_score = NA_real_,
      centroid_dist = NA_real_,
      concentration = NA_real_,
      peak = NA_real_
    ))
  }
  
  # hotspot uses positive mass only (where B-A indicates "worse" for example)
  w <- pmax(d$delta, 0)
  if (sum(w) <= 0) {
    return(list(
      hotspot_score = NA_real_,
      centroid_dist = 0,
      concentration = NA_real_,
      peak = max(d$delta, na.rm=TRUE)
    ))
  }
  
  xcm <- sum(d$Xb * w) / sum(w)
  ycm <- sum(d$Yb * w) / sum(w)
  
  dist <- sqrt((xcm - cx)^2 + (ycm - cy)^2)
  
  # concentration: weighted 2nd moment around centroid (smaller => more concentrated hotspot)
  r2 <- (d$Xb - xcm)^2 + (d$Yb - ycm)^2
  m2 <- sum(w * r2) / sum(w)
  conc <- 1 / sqrt(m2 + 1e-9)  # bigger => more concentrated
  
  peak <- max(d$delta, na.rm=TRUE)
  
  # composite score (tunable)
  # - dist: hotspot is meaningful if off-center
  # - conc: localized hotspot is more "hotspot-like"
  # - peak: strong signal
  score <- dist * conc * abs(peak)
  
  list(
    hotspot_score = score,
    centroid_dist = dist,
    concentration = conc,
    peak = peak
  )
}

# -----------------------------
# compute ranking over all MSRs
# -----------------------------
# -----------------------------
# compute ranking over all MSRs (+ progress bar)
# -----------------------------
rank_list <- vector("list", length(msr_names))

pb <- txtProgressBar(min = 0, max = length(msr_names), style = 3)
for (i in seq_along(msr_names)) {
  msr <- msr_names[i]
  
  # group tile mean for this MSR
  tile <- dt[, .(mean_score = mean(get(msr), na.rm=TRUE)),
             by=.(GROUP, Xb, Yb)]
  
  # complete grid for A and B
  Amap <- merge(grid, tile[GROUP=="A", .(Xb,Yb, A=mean_score)], by=c("Xb","Yb"), all.x=TRUE)
  Bmap <- merge(grid, tile[GROUP=="B", .(Xb,Yb, B=mean_score)], by=c("Xb","Yb"), all.x=TRUE)
  
  Dmap <- merge(Amap, Bmap, by=c("Xb","Yb"), all=TRUE)
  Dmap[, delta := B - A]
  
  met <- hotspot_metrics(Dmap, cx=CX, cy=CY)
  
  rank_list[[i]] <- data.table(
    MSR = msr,
    hotspot_score = met$hotspot_score,
    centroid_dist = met$centroid_dist,
    concentration = met$concentration,
    peak = met$peak
  )
  
  setTxtProgressBar(pb, i)
}
close(pb)


rank_dt <- rbindlist(rank_list)
rank_dt <- rank_dt[!is.na(hotspot_score)]
setorder(rank_dt, -hotspot_score)

# show top 10 table
top10 <- head(rank_dt, 10)
print(top10)

# -----------------------------
# gallery plots: for each MSR in top10, print A / B / Delta maps
# -----------------------------
plot_one_msr <- function(msr) {
  tile <- dt[, .(mean_score = mean(get(msr), na.rm=TRUE)),
             by=.(GROUP, Xb, Yb)]
  
  Amap <- merge(grid, tile[GROUP=="A", .(Xb,Yb, val=mean_score)], by=c("Xb","Yb"), all.x=TRUE)
  Bmap <- merge(grid, tile[GROUP=="B", .(Xb,Yb, val=mean_score)], by=c("Xb","Yb"), all.x=TRUE)
  
  Dmap <- merge(Amap[,.(Xb,Yb,A=val)], Bmap[,.(Xb,Yb,B=val)], by=c("Xb","Yb"), all=TRUE)
  Dmap[, delta := B - A]
  
  pA <- ggplot(Amap, aes(Xb,Yb,fill=val)) +
    geom_tile(width=BIN,height=BIN) +
    geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.35) +
    coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
    theme_minimal(base_size=11) +
    labs(title=paste0(msr, " | GROUP A"), fill="mean")
  
  pB <- ggplot(Bmap, aes(Xb,Yb,fill=val)) +
    geom_tile(width=BIN,height=BIN) +
    geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.35) +
    coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
    theme_minimal(base_size=11) +
    labs(title=paste0(msr, " | GROUP B"), fill="mean")
  
  pD <- ggplot(Dmap, aes(Xb,Yb,fill=delta)) +
    geom_tile(width=BIN,height=BIN) +
    geom_path(data=circle, aes(x=x,y=y), inherit.aes=FALSE, linewidth=0.35) +
    coord_equal(xlim=c(CX-R,CX+R), ylim=c(CY-R,CY+R)) +
    theme_minimal(base_size=11) +
    labs(title=paste0(msr, " | Δ (B-A)"), fill="delta")
  
  print(pA); print(pB); print(pD)
}

pb2 <- txtProgressBar(min = 0, max = nrow(top10), style = 3)
for (j in seq_len(nrow(top10))) {
  plot_one_msr(top10$MSR[j])
  setTxtProgressBar(pb2, j)
}
close(pb2)

