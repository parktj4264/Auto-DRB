# src/01_load.R ------------------------------------------------------------


# 1. 파일 경로 설정 및 체크 -------------------------------------------------
cat("Reading data...\n")

path_raw  <- here::here("data", FILE_NAME_RAW)
path_meta <- here::here("data", FILE_NAME_META)

if (!file.exists(path_raw))  stop(paste("Raw Data 없음:", path_raw))
if (!file.exists(path_meta)) stop(paste("Meta Data 없음:", path_meta))


# 2. 데이터 로드 (Fast Read) ------------------------------------------------
dt_raw  <- fread(path_raw,  header = TRUE, nThread = getDTthreads())
dt_meta <- fread(path_meta, header = TRUE, nThread = getDTthreads())

cat(sprintf("Raw Loaded : %s rows, %s cols\n", format(nrow(dt_raw), big.mark=","), ncol(dt_raw)))
cat(sprintf("Meta Loaded: %s rows, %s cols\n", format(nrow(dt_meta), big.mark=","), ncol(dt_meta)))


# 3. 필수 컬럼 유효성 체크 ----------------------------------------------------
# raw 필수: ROOTID, PARTID
req_raw  <- c(ROOTID_COL, PARTID_COL)
miss_raw <- setdiff(req_raw, names(dt_raw))
if (length(miss_raw) > 0) {
  stop(sprintf("raw.csv에 필수 컬럼 없음: %s", paste(miss_raw, collapse = ", ")))
}

# meta 필수: ROOTID, GROUP
req_meta  <- c(ROOTID_COL, GROUP_COL)
miss_meta <- setdiff(req_meta, names(dt_meta))
if (length(miss_meta) > 0) {
  stop(sprintf("ROOTID.csv에 필수 컬럼 없음: %s", paste(miss_meta, collapse = ", ")))
}


# 4. MSR 컬럼 식별 (PARTID 이후 전부) ---------------------------------------
cat(">> Identifying MSR columns based on 'PARTID'...\n")

idx_partid <- which(names(dt_raw) == PARTID_COL)
if (length(idx_partid) != 1) stop("raw.csv에서 PARTID 컬럼이 없거나(또는 중복)합니다.")

idx_start <- idx_partid + 1
if (idx_start > ncol(dt_raw)) stop("PARTID가 마지막 컬럼입니다. PARTID 뒤에 MSR 컬럼이 없습니다.")

MSR_COLS <- names(dt_raw)[idx_start:ncol(dt_raw)]
if (length(MSR_COLS) == 0) stop("MSR 컬럼이 0개입니다. (PARTID 뒤에 데이터 없음)")

cat(sprintf("Target MSR Count: %d\n", length(MSR_COLS)))


# 5. meta GROUP 값 체크 (REF/TARGET 존재 여부) ------------------------------
# (원칙: run.R의 REF/TARGET를 존중하되, 불일치하면 ROOTID.csv 기준으로 자동 지정)
dt_meta[, (GROUP_COL) := as.character(get(GROUP_COL))]

u_groups <- sort(unique(dt_meta[[GROUP_COL]]))
cat(sprintf(">> Meta GROUP levels: %s\n", paste(u_groups, collapse = ", ")))

if (length(u_groups) < 2) {
  stop(sprintf("Meta의 GROUP 유니크 값이 2개 미만입니다. (현재: %s)  REF/TARGET 비교 불가",
               paste(u_groups, collapse=", ")))
}

# run.R 설정 라벨이 meta에 실제로 존재하는지 확인
ok_ref    <- (GROUP_REF_LABEL    %in% u_groups)
ok_target <- (GROUP_TARGET_LABEL %in% u_groups)

if (!ok_ref || !ok_target) {
  cat("WARNING: run.R의 REF/TARGET 라벨이 meta와 불일치합니다. ROOTID.csv 기준으로 REF/TARGET를 자동 지정합니다.\n")
  cat(sprintf(" - run.R REF/TARGET : %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
  
  # 자동 지정 규칙: 정렬된 유니크 값의 1번=REF, 2번=TARGET
  GROUP_REF_LABEL    <- u_groups[1]
  GROUP_TARGET_LABEL <- u_groups[2]
  
  cat(sprintf(" - auto  REF/TARGET : %s / %s\n", GROUP_REF_LABEL, GROUP_TARGET_LABEL))
}



# 6. 데이터 병합 (Merge): raw에 GROUP 붙이기 --------------------------------
cat(">> Merging raw + meta by ROOTID...\n")

# merge 결과 컬럼 순서가 바뀔 수 있으니 MSR_COLS는 위에서 raw 기준으로 이미 확정!
dt <- merge(dt_raw, dt_meta[, c(ROOTID_COL, GROUP_COL), with=FALSE],
            by = ROOTID_COL, all.x = TRUE)

# GROUP 누락 체크 (stop 대신 warning + 제외)
n_na_group <- sum(is.na(dt[[GROUP_COL]]))
if (n_na_group > 0) {
  cat(sprintf("WARNING: Merge 후 GROUP이 NA인 row가 %s개 있습니다. 해당 row는 제외하고 진행합니다.\n",
              format(n_na_group, big.mark=",")))
  
  # 어떤 ROOTID가 누락됐는지 샘플 표시
  miss_ids <- unique(dt[is.na(get(GROUP_COL)), get(ROOTID_COL)])
  cat(sprintf(" - 누락 ROOTID 예시: %s\n",
              paste(head(miss_ids, 5), collapse=", ")))
  
  # 제외
  dt <- dt[!is.na(get(GROUP_COL))]
}


# GROUP 타입 통일
dt[, (GROUP_COL) := as.character(get(GROUP_COL))]


# 7. MSR 컬럼 타입 정리 (numeric 강제) ---------------------------------------
# 실무 데이터에서 MSR에 문자/공백 섞이는 사고 방지
cat(">> Coercing MSR columns to numeric (safe)...\n")

for (nm in MSR_COLS) {
  if (!is.numeric(dt[[nm]])) {
    suppressWarnings({
      dt[, (nm) := as.numeric(get(nm))]
    })
  }
}


# 8. 완료 로그 --------------------------------------------------------------
cat("Done (Merged)\n")
cat(sprintf("Final Data: %s rows, %s cols\n", format(nrow(dt), big.mark=","), ncol(dt)))
cat(sprintf("Group column: %s (REF=%s, TARGET=%s)\n", GROUP_COL, GROUP_REF_LABEL, GROUP_TARGET_LABEL))
# cat(sprintf("MSR range: %s ... %s\n", MSR_COLS[1], MSR_COLS[length(MSR_COLS)]))
