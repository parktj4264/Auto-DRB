# 04_save.R
# src/04_save.R ------------------------------------------------------------

cat("Saving outputs...\n")

# 1) 결과 객체 체크 ----------------------------------------------------------
if (!exists("results_dt")) {
  cat("WARNING: results_dt가 없습니다. 저장을 스킵합니다.\n")
} else {
  
  # 2) output dir 준비 ------------------------------------------------------
  out_dir <- here::here(OUT_DIR)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf(">> Created output dir: %s\n", out_dir))
  }
  
  # 3) 파일 저장 ------------------------------------------------------------
  out_path <- file.path(out_dir, OUT_FILE_RESULTS)
  
  # overwrite 안내 (실무에서 중요)
  if (file.exists(out_path)) {
    cat(sprintf("WARNING: existing file will be overwritten: %s\n", out_path))
  }
  
  # [NEW] flag는 CSV에서 1/0으로 저장 (실무 편의)
  # [FIX] flag는 CSV에서 1/0으로 저장 (TRUE=1, FALSE=0, NA=0)
  flag_cols <- intersect(c("sigma_up","sigma_down","wilcox_flag","ks_flag"), names(results_dt))
  if (length(flag_cols) > 0) {
    for (fc in flag_cols) {
      results_dt[, (fc) := {
        v <- as.integer(get(fc))  # TRUE=1, FALSE=0, NA=NA
        v[is.na(v)] <- 0L         # NA는 0 처리
        v
      }]
    }
  }
  
  
  # fwrite (빠르고 안전)
  tryCatch({
    data.table::fwrite(results_dt, out_path)
    cat(sprintf(">> Saved: %s\n", out_path))
  }, error = function(e) {
    cat("ERROR: fwrite failed.\n")
    cat(sprintf(" - message: %s\n", e$message))
  })
  
  # 4) 간단 요약 로그 -------------------------------------------------------
  cat(sprintf(">> rows: %s | cols: %d\n",
              format(nrow(results_dt), big.mark=","), ncol(results_dt)))
  
  # flag 요약 (0/1이면 sum이 곧 count)
  if (all(c("sigma_up","sigma_down","wilcox_flag","ks_flag") %in% names(results_dt))) {
    cat(sprintf(">> flags: sigma_up=%d, sigma_down=%d, wilcox=%d, ks=%d\n",
                sum(results_dt$sigma_up, na.rm=TRUE),
                sum(results_dt$sigma_down, na.rm=TRUE),
                sum(results_dt$wilcox_flag, na.rm=TRUE),
                sum(results_dt$ks_flag, na.rm=TRUE)))
  }
}

cat("Done (Save)\n")

