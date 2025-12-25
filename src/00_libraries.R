# 패키지 자동 설치 및 로딩 함수 (무적 버전) -------------------------------------------------------
library_load <- function(packages) {
  
  # [핵심] 옵션 강제 설정: 질문 금지, 바이너리 강제, 주소 지정
  options(repos = c(CRAN = "https://cran.rstudio.com/"))  # 다운로드 주소 고정
  options(pkgType = "win.binary")                         # 윈도우용 완제품만 (컴파일 X)
  options(install.packages.check.source = "no")           # 소스 확인 안함
  options(install.packages.compile.from.source = "never") # 컴파일 절대 안함 (에러 방지)
  
  # 색상 정의 (콘솔 로그 가독성 UP)
  green  <- function(x) paste0("\033[32m", x, "\033[0m")
  yellow <- function(x) paste0("\033[33m", x, "\033[0m")
  blue   <- function(x) paste0("\033[34m", x, "\033[0m")
  red    <- function(x) paste0("\033[31m", x, "\033[0m")
  gray   <- function(x) paste0("\033[90m", x, "\033[0m")
  
  total <- length(packages)
  
  for (i in seq_along(packages)) {
    package <- packages[i]
    message(gray(strrep("-", 50)))
    message(gray(paste0("Package [", i, "/", total, "]")))
    
    if (!requireNamespace(package, quietly = TRUE)) {
      message(yellow(paste("Installing:", package, "(Binary Only)")))
      
      tryCatch({
        # 여기서 type="binary" 한번 더 명시해서 확인사살
        install.packages(package, type = "binary", quiet = TRUE)
      }, error = function(e) {
        message(red(paste("Install failed:", package)))
        message(red(paste("Error:", e$message)))
      })
      
    } else {
      message(green(paste("Already installed:", package)))
    }
    
    message(blue(paste("Loading:", package)))
    suppressPackageStartupMessages(
      library(package, character.only = TRUE)
    )
  }
  
  message(gray(strrep("-", 50)))
  message(green("All requested packages processed."))
}

# 사용할 패키지 목록 --------------------------------------------------------------
cat("library를 불러옵니다...\n")

library_load(
  c("data.table", "here", "stringr", "lubridate", "purrr", "stats")
)