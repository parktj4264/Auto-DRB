library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
dt_raw   <- fread("data/raw.csv")
dt_meta  <- fread("data/ROOTID.csv")
dt_truth <- fread("data/truth_msr.csv")
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

# 랜덤 MSR 하나 뽑기 (hotspot 우선, 없으면 mean)
pick_one <- function(effect) {
  x <- dt_truth[effect_type == effect, MSR]
  if (length(x) == 0) return(NA_character_)
  sample(x, 1)
}
# set.seed(123)
msr <- pick_one("hotspot")
if (is.na(msr)) msr <- pick_one("mean")
stopifnot(!is.na(msr))
cat("Using MSR =", msr, "\n")

# ------------------------------------------------------------
# Helpers: Gaussian kernel + conv2(same)
# ------------------------------------------------------------
make_gaussian_kernel <- function(sigma = 1.0) {
  if (!is.finite(sigma) || sigma <= 0) return(matrix(1, 1, 1))
  k_size <- ceiling(3 * sigma) * 2 + 1
  center <- (k_size + 1) / 2
  ii <- seq_len(k_size)
  d2 <- (ii - center)^2
  kernel <- exp(-(outer(d2, d2, "+")) / (2 * sigma^2))
  kernel / sum(kernel)
}

conv2_same <- function(mat, kernel) {
  nr <- nrow(mat); nc <- ncol(mat)
  kr <- nrow(kernel); kc <- ncol(kernel)
  khr <- floor(kr / 2); khc <- floor(kc / 2)
  
  padded <- matrix(0, nr + 2*khr, nc + 2*khc)
  padded[(khr+1):(nr+khr), (khc+1):(nc+khc)] <- mat
  
  out <- matrix(0, nr, nc)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      sub_m <- padded[i:(i+2*khr), j:(j+2*khc)]
      out[i, j] <- sum(sub_m * kernel)
    }
  }
  out
}

# ------------------------------------------------------------
# Connected components (4-neighborhood) + size filter
#   bin: logical matrix -> keep logical matrix
# ------------------------------------------------------------
filter_small_components4 <- function(bin, min_size = 10L) {
  nr <- nrow(bin); nc <- ncol(bin)
  visited <- matrix(FALSE, nr, nc)
  keep <- matrix(FALSE, nr, nc)
  
  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      if (!bin[r, c] || visited[r, c]) next
      
      q_r <- c(r); q_c <- c(c)
      visited[r, c] <- TRUE
      
      comp_r <- integer(0); comp_c <- integer(0)
      head <- 1L
      
      while (head <= length(q_r)) {
        rr <- q_r[head]; cc <- q_c[head]
        head <- head + 1L
        
        comp_r <- c(comp_r, rr)
        comp_c <- c(comp_c, cc)
        
        if (rr > 1L && bin[rr-1L, cc] && !visited[rr-1L, cc]) {
          visited[rr-1L, cc] <- TRUE; q_r <- c(q_r, rr-1L); q_c <- c(q_c, cc)
        }
        if (rr < nr && bin[rr+1L, cc] && !visited[rr+1L, cc]) {
          visited[rr+1L, cc] <- TRUE; q_r <- c(q_r, rr+1L); q_c <- c(q_c, cc)
        }
        if (cc > 1L && bin[rr, cc-1L] && !visited[rr, cc-1L]) {
          visited[rr, cc-1L] <- TRUE; q_r <- c(q_r, rr); q_c <- c(q_c, cc-1L)
        }
        if (cc < nc && bin[rr, cc+1L] && !visited[rr, cc+1L]) {
          visited[rr, cc+1L] <- TRUE; q_r <- c(q_r, rr); q_c <- c(q_c, cc+1L)
        }
      }
      
      if (length(comp_r) >= min_size) keep[cbind(comp_r, comp_c)] <- TRUE
    }
  }
  keep
}

# ------------------------------------------------------------
# Stage builder: raw -> z-clip -> smooth -> soft-binarize -> CC filter
#   map_dt: data.table(x,y,v)
# ------------------------------------------------------------
calc_stages_soft_cc <- function(map_dt,
                                sigma_thresh = 3.0,
                                smooth_sigma = 1.5,
                                tau_q = 0.90,
                                min_size = 12L,
                                tiny = 1e-12) {
  dt0 <- as.data.table(map_dt)[, .(x, y, v)]
  dt0 <- dt0[is.finite(x) & is.finite(y)]
  dt0 <- dt0[, .(v = mean(v, na.rm=TRUE)), by=.(x,y)]
  
  xs <- sort(unique(dt0$x))
  ys <- sort(unique(dt0$y))
  x_to_c <- match(dt0$x, xs)
  y_to_r <- match(dt0$y, ys)
  
  # raw matrix (for display)
  mat_raw <- matrix(0, nrow=length(ys), ncol=length(xs))
  mat_raw[cbind(y_to_r, x_to_c)] <- dt0$v
  
  # Z
  vals <- dt0$v
  mu <- mean(vals, na.rm=TRUE)
  sdv <- sd(vals, na.rm=TRUE)
  
  mat_w0 <- matrix(0, nrow=length(ys), ncol=length(xs))
  if (is.finite(sdv) && sdv >= tiny) {
    z <- (vals - mu) / sdv
    z[!is.finite(z)] <- 0
    
    # Z-filter/clip (추천: clip-relu)
    # abs(z)에서 thresh만큼 잘라내고 남은 것만 에너지로
    w0 <- pmax(abs(z) - sigma_thresh, 0)
    mat_w0[cbind(y_to_r, x_to_c)] <- w0
  }
  
  # Smooth
  kernel <- make_gaussian_kernel(smooth_sigma)
  mat_sm <- conv2_same(mat_w0, kernel)
  
  # Soft-binarize (등고선 느낌): 상위 q 분위수 기준으로 아래는 0으로 뭉개기
  tau <- as.numeric(quantile(mat_sm, probs=tau_q, na.rm=TRUE, names=FALSE, type=7))
  if (!is.finite(tau)) tau <- 0
  mat_soft <- pmax(mat_sm - tau, 0)
  
  # Connected component filter (고립 제거): soft>0의 덩어리 중 작은 덩어리 삭제
  bin <- mat_soft > 0
  keep <- filter_small_components4(bin, min_size = min_size)
  mat_final <- mat_soft * keep
  
  # Build long dt for ggplot (per-point)
  dt_out <- data.table(
    x = dt0$x, y = dt0$y,
    raw   = mat_raw[cbind(y_to_r, x_to_c)],
    w0    = mat_w0[cbind(y_to_r, x_to_c)],
    sm    = mat_sm[cbind(y_to_r, x_to_c)],
    soft  = mat_soft[cbind(y_to_r, x_to_c)],
    final = mat_final[cbind(y_to_r, x_to_c)]
  )
  
  out <- melt(dt_out, id.vars=c("x","y"),
              measure.vars=c("raw","w0","sm","soft","final"),
              variable.name="stage", value.name="val")
  out[, val := fifelse(!is.finite(val), 0, val)]
  
  # stage-wise scale to [0,1] for fair visual comparison
  out[, val_scaled := {
    mn <- min(val, na.rm=TRUE); mx <- max(val, na.rm=TRUE)
    if (!is.finite(mx - mn) || (mx - mn) < tiny) rep(0, .N) else (val - mn) / (mx - mn)
  }, by=stage]
  
  out[, stage := factor(stage,
                        levels=c("raw","w0","sm","soft","final"),
                        labels=c("raw (GROUP mean map)",
                                 "Z-filter/clip: pmax(|z|-t,0)",
                                 "Smooth(w0)",
                                 sprintf("soft-binarize: pmax(sm-q%.2f,0)", tau_q),
                                 sprintf("CC filter: keep size >= %d", min_size)))]
  out[]
}

# ------------------------------------------------------------
# Make GROUP mean map (raw)
# ------------------------------------------------------------
dt_map <- dt[, .(v = mean(get(msr), na.rm=TRUE)), by=.(GROUP, X, Y)]
setnames(dt_map, c("X","Y"), c("x","y"))

# Parameters (튜닝 포인트)
sigma_thresh <- 1.0
smooth_sigma <- 1.5
tau_q        <- 0.90
min_size     <- 12L

stages <- dt_map[, calc_stages_soft_cc(.SD,
                                       sigma_thresh=sigma_thresh,
                                       smooth_sigma=smooth_sigma,
                                       tau_q=tau_q,
                                       min_size=min_size),
                 by=.(GROUP)]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(stages, aes(x=x, y=y, fill=val_scaled)) +
  geom_tile() +
  coord_fixed() +
  facet_grid(stage ~ GROUP) +
  scale_fill_viridis_c(option="C", direction=1) +
  labs(
    title = sprintf("AB what-smoothing-sees (soft-binarize + CC) | %s | zt=%.1f, smooth=%.1f, q=%.2f, min_size=%d",
                    msr, sigma_thresh, smooth_sigma, tau_q, min_size),
    x="X", y="Y", fill="scaled"
  ) +
  theme_minimal(base_size=12)

print(p)
