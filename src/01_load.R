# src/01_load.R
# Written by: Tae-jun Park
# Description: Raw Data 및 Meta Data 로드, 병합, MSR 컬럼 자동 식별

# 1. 파일 경로 설정 및 체크 --------------------------------------------------------
cat(green(">> [Step 1] 파일 경로 확인 중...\n"))

path_raw  <- here::here("data", FILE_NAME_RAW)
path_meta <- here::here("data", FILE_NAME_META)

# 파일 존재 여부 확인 (Samsung Style: 데이터 없으면 바로 Stop)
if (!file.exists(path_raw))  stop(red(sprintf("Raw Data가 없습니다: %s", path_raw)))
if (!file.exists(path_meta)) stop(red(sprintf("Meta Data가 없습니다: %s", path_meta)))

cat(sprintf("   - Raw Path : %s\n", basename(path_raw)))
cat(sprintf("   - Meta Path: %s\n", basename(path_meta)))


# 2. 데이터 로드 (Fast Read by data.table) -------------------------------------
cat(green("\n>> [Step 2] 데이터 고속 로드 (fread)...\n"))

# 스레드 수 명시 (형 워크스테이션 성능 다 뽑아먹자)
dt_raw  <- fread(path_raw,  header = TRUE, nThread = getDTthreads())
dt_meta <- fread(path_meta, header = TRUE, nThread = getDTthreads())

cat(sprintf("   - Raw Loaded : %s rows, %s cols\n", format(nrow(dt_raw), big.mark=","), ncol(dt_raw)))
cat(sprintf("   - Meta Loaded: %s rows, %s cols\n", format(nrow(dt_meta), big.mark=","), ncol(dt_meta)))


# 3. 유효성 체크 및 그룹 변수 정리 -------------------------------------------------
cat(green("\n>> [Step 3] 그룹 변수 유효성 검증...\n"))

if (length(ADD_GROUPS) > 0) {
  # 메타 데이터에 실제로 존재하는 컬럼만 교집합(intersect)으로 살림
  valid_add_groups <- intersect(ADD_GROUPS, names(dt_meta))
  missing_cols     <- setdiff(ADD_GROUPS, names(dt_meta))
  
  # 없는 컬럼은 경고(warning)만 주고 쿨하게 패스
  if (length(missing_cols) > 0) {
    message(yellow(sprintf("   ⚠️ 경고: 요청한 그룹 중 '%s' 컬럼이 Meta에 없습니다. -> 제외합니다.", 
                           paste(missing_cols, collapse=", "))))
  }
  
  # ADD_GROUPS 최신화
  ADD_GROUPS <- valid_add_groups
}

# 최종 그룹핑 컬럼 확정 (ROOTID + GROUP + 살아남은 ADD_GROUPS)
GROUP_COLS <- c("ROOTID", "GROUP", ADD_GROUPS)
cat(sprintf("   - Final Group Cols: %s\n", paste(GROUP_COLS, collapse = ", ")))


# 4. 데이터 병합 (Left Join) -----------------------------------------------------
cat(green("\n>> [Step 4] Raw + Meta 병합 (Left Join)...\n"))

# Raw 데이터 기준 병합 (all.x = TRUE 필수, 메타 없어도 측정값은 살려야지)
dt <- merge(dt_raw, dt_meta, by = "ROOTID", all.x = TRUE)

# 메모리 정리 (Tiny Data 철학: 다 쓴 객체는 바로바로 삭제)
rm(dt_raw, dt_meta); gc() 


# 5. MSR 컬럼 식별 (Dynamic Logic) ---------------------------------------------
cat(green("\n>> [Step 5] 'PARTID' 기준 MSR 데이터 컬럼 식별...\n"))

if ("PARTID" %in% names(dt)) {
  
  # PARTID 컬럼의 인덱스 찾기
  idx_partid <- which(names(dt) == "PARTID")
  
  # PARTID가 마지막 컬럼이면 뒤에 데이터가 없는 거니까 에러 처리
  if (idx_partid == ncol(dt)) {
    stop(red("'PARTID'가 마지막 컬럼입니다. 뒤에 MSR 데이터가 존재하지 않습니다."))
  }
  
  # PARTID 바로 뒤부터 끝까지가 MSR 데이터
  MSR_COLS <- names(dt)[(idx_partid + 1):ncol(dt)]
  
  # [선택] 숫자형(Numeric) 체크 로직 (필요하면 주석 해제)
  # MSR 데이터에 'Fail', 'Error' 같은 문자열 섞여있으면 여기서 걸러야 함
  # numeric_cols <- sapply(dt[, ..MSR_COLS], is.numeric)
  # if (!all(numeric_cols)) {
  #   message(yellow(sprintf("주의: MSR 컬럼 중 %d개가 숫자가 아닙니다 (제외됨)", sum(!numeric_cols))))
  #   MSR_COLS <- MSR_COLS[numeric_cols]
  # }
  
  if (length(MSR_COLS) == 0) stop(red("유효한 MSR 컬럼이 0개입니다."))
  
} else {
  # PARTID 없으면 분석 불가
  stop(red("데이터에 'PARTID' 컬럼이 없습니다. (Structure Error)"))
}


# 완료 리포트 ------------------------------------------------------------------
cat(gray(strrep("-", 60)), "\n")
cat(green("[Success] Data Load & Merge Completed!\n"))
cat(sprintf("   - Final Table Shape : %s rows x %s cols\n", 
            format(nrow(dt), big.mark=","), ncol(dt)))
cat(sprintf("   - Target MSR Count  : %d 개\n", length(MSR_COLS)))
cat(sprintf("   - Memory Usage      : %s\n", format(object.size(dt), units="auto")))
cat(gray(strrep("-", 60)), "\n")




