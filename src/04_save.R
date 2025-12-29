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
  
  # [OPTIONAL] direction을 코드로도 남기고 싶으면 추가 (Up=1, Stable=0, Down=-1)
  # - 기본은 direction 문자열 유지 (엔지니어 가독성)
  if ("direction" %in% names(results_dt) && !("dir_code" %in% names(results_dt))) {
    results_dt[, dir_code := fifelse(direction == "Up", 1L,
                                     fifelse(direction == "Down", -1L, 0L))]
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
