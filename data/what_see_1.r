library(data.table)
library(ggplot2)

dt_raw   <- fread("data/raw.csv")
dt_meta  <- fread("data/ROOTID.csv")
dt_truth <- fread("data/truth_msr.csv")
dt <- merge(dt_raw, dt_meta, by="ROOTID", all.x=TRUE)

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

calc_smoothing_stages_with_raw <- function(map_dt,
                                           sigma_thresh = 3.0,
                                           smooth_sigma = 1.5,
                                           iso_alpha = 2.0,
                                           tiny = 1e-12) {
  dt0 <- as.data.table(map_dt)[, .(x, y, v)]
  dt0 <- dt0[is.finite(x) & is.finite(y)]
  dt0 <- dt0[, .(v = mean(v, na.rm=TRUE)), by=.(x,y)]
  
  vals <- dt0$v
  mu <- mean(vals, na.rm=TRUE)
  sdv <- sd(vals, na.rm=TRUE)
  if (!is.finite(sdv) || sdv < tiny) {
    dt0[, `:=`(raw=v, w0=0, sm=0, d=0, final=0)]
  } else {
    z <- (vals - mu) / sdv
    z[!is.finite(z)] <- 0
    
    if (!is.finite(sigma_thresh) || sigma_thresh < 0) sigma_thresh <- 0
    w0 <- ifelse(abs(z) > sigma_thresh, abs(z), 0)
    
    xs <- sort(unique(dt0$x))
    ys <- sort(unique(dt0$y))
    x_to_c <- match(dt0$x, xs)
    y_to_r <- match(dt0$y, ys)
    
    mat_w0 <- matrix(0, nrow=length(ys), ncol=length(xs))
    mat_w0[cbind(y_to_r, x_to_c)] <- w0
    mat_mask <- (mat_w0 > 0) * 1.0
    
    kernel <- make_gaussian_kernel(smooth_sigma)
    sm_mat <- conv2_same(mat_w0, kernel)
    d_mat  <- conv2_same(mat_mask, kernel)
    final_mat <- sm_mat * (d_mat ^ iso_alpha)
    
    dt0[, raw := v]
    dt0[, w0 := w0]
    dt0[, sm := as.numeric(sm_mat[cbind(y_to_r, x_to_c)])]
    dt0[, d  := as.numeric(d_mat[cbind(y_to_r, x_to_c)])]
    dt0[, final := as.numeric(final_mat[cbind(y_to_r, x_to_c)])]
  }
  
  out <- melt(
    dt0,
    id.vars=c("x","y"),
    measure.vars=c("raw","w0","sm","d","final"),
    variable.name="stage",
    value.name="val"
  )
  
  out[, val := fifelse(!is.finite(val), 0, val)]
  out[, val_scaled := {
    mn <- min(val, na.rm=TRUE)
    mx <- max(val, na.rm=TRUE)
    if (!is.finite(mx - mn) || (mx - mn) < tiny) rep(0, .N) else (val - mn) / (mx - mn)
  }, by=stage]
  
  out[]
}

# GROUP mean map for chosen MSR
dt_map <- dt[, .(v = mean(get(msr), na.rm=TRUE)), by=.(GROUP, X, Y)]
setnames(dt_map, c("X","Y"), c("x","y"))

sigma_thresh <- 1.0
smooth_sigma <- 1.5
iso_alpha    <- 1.5

stages <- dt_map[, calc_smoothing_stages_with_raw(.SD,
                                                  sigma_thresh=sigma_thresh,
                                                  smooth_sigma=smooth_sigma,
                                                  iso_alpha=iso_alpha),
                 by=.(GROUP)]

stages[, stage := factor(stage,
                         levels=c("raw","w0","sm","d","final"),
                         labels=c("raw (GROUP mean map)",
                                  "w0 (Z-filter energy)",
                                  "Smooth(w0)",
                                  "d = Smooth(mask)",
                                  "final = Smooth(w0) * d^alpha"))]

p <- ggplot(stages, aes(x=x, y=y, fill=val_scaled)) +
  geom_tile() +
  coord_fixed() +
  facet_grid(stage ~ GROUP) +
  scale_fill_viridis_c(option="C", direction=1) +
  labs(
    title = sprintf("What smoothing sees (with raw) | %s | sigma=%.1f, smooth=%.1f, alpha=%.1f",
                    msr, sigma_thresh, smooth_sigma, iso_alpha),
    x="X", y="Y", fill="scaled"
  ) +
  theme_minimal(base_size=12)

print(p)
