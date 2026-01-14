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

  # overwrite 안내
  if (file.exists(out_path)) {
    rel_out_warn <- file.path(OUT_DIR, OUT_FILE_RESULTS)
    cat(sprintf("WARNING: existing file will be overwritten: %s\n", rel_out_warn))
  }

  # [OPTIONAL] direction code
  if ("direction" %in% names(results_dt) && !("dir_code" %in% names(results_dt))) {
    results_dt[, dir_code := fifelse(
      direction == "Up", 1L,
      fifelse(direction == "Down", -1L, 0L)
    )]
  }

  # (A) Main Save (Overwrite)
  tryCatch(
    {
      data.table::fwrite(results_dt, out_path)
      # Relative path for display
      rel_out <- file.path(OUT_DIR, OUT_FILE_RESULTS)
      cat(sprintf(">> Saved (Main): %s\n", rel_out))
    },
    error = function(e) {
      cat("ERROR: Main fwrite failed.\n")
      cat(sprintf(" - message: %s\n", e$message))
    }
  )

  # (B) History Save (Timestamped)
  tryCatch(
    {
      # timestamp from log filename or current time
      # (일관성을 위해 현재 시간 다시 찍음)
      ts_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
      hist_dir <- file.path(out_dir, ts_str)

      if (!dir.exists(hist_dir)) dir.create(hist_dir)

      # Filename also timestamped
      hist_file <- sprintf("results_%s.csv", ts_str)
      hist_path <- file.path(hist_dir, hist_file)

      data.table::fwrite(results_dt, hist_path)

      # Relative path for display
      rel_hist <- file.path(OUT_DIR, ts_str, hist_file)
      cat(sprintf(">> Saved (History): %s\n", rel_hist))
    },
    error = function(e) {
      cat("WARNING: History save failed (Main save is safe).\n")
      cat(sprintf(" - message: %s\n", e$message))
    }
  )

  # 4) 간단 요약 로그 -------------------------------------------------------
  cat(sprintf(
    ">> rows: %s | cols: %d\n",
    format(nrow(results_dt), big.mark = ","), ncol(results_dt)
  ))

  # direction 요약
  if ("direction" %in% names(results_dt)) {
    cat(">> direction summary:\n")
    print(table(results_dt$direction, useNA = "ifany"))
  }

  # # score 요약 (엔지니어가 바로 감 잡음)
  # score_cols <- intersect(c("sigma_score", "cliffs_delta", "ws_spatial"), names(results_dt))
  # if (length(score_cols) > 0) {
  #   cat(">> score summary (quantiles):\n")
  #   for (sc in score_cols) {
  #     v <- results_dt[[sc]]
  #     v <- v[is.finite(v)]
  #     if (length(v) == 0) next
  #     qs <- stats::quantile(v, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
  #     cat(sprintf(" - %s: ", sc))
  #     print(qs)
  #   }
  # }

  # TOP10 예시 출력 (abs(sigma_score) 기준으로 이미 정렬돼있다는 가정)
  # if ("sigma_score" %in% names(results_dt)) {
  #   cat(">> TOP 10 by |sigma_score| (preview):\n")
  #   cols_show <- intersect(c("msr","direction","sigma_score","cliffs_delta","ws_spatial",
  #                            "mean_ref","sd_ref","mean_target","sd_target"), names(results_dt))
  #   print(head(results_dt[, ..cols_show], 10))
  # }
}

cat("Done (Save)\n")
