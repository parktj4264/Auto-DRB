# DRB rev1.0 최종 목표

-   A공정(REF) → B공정(TARGET)으로 바뀔 때, 흔들리는 MSR을 자동으로 뽑아준다. 이 결과는 **판정(자동 결론)**이 아니라, 엔지니어가 눈으로 맵/분포를 확인할 후보 리스트를 만들어주는 “리스트업 필터”다.
-   DRB rev1.0은 REF 공정을 기준으로 TARGET 공정의 레벨·분포·공간 성향 변화를 탐지하는 리뷰 보조 도구이다.

## 데이터 전제 (실데이터 기준)

-   `raw.csv` : 칩 레벨 데이터
    -   키: `ROOTID`
    -   필수: `PARTID`, `Radius`
    -   MSR 컬럼: `PARTID` 바로 다음 열부터 끝까지 전부
-   `ROOTID.csv` : 웨이퍼 레벨 메타
    -   `ROOTID`, `GROUP` (A/B)

> **주의**: 실무의 `Radius`는 이름만 Radius고, 실제 의미는 좌/우 signed position (-150 \~ 150, 실수)로 가정한다.

## DRB rev1.0 분석 렌즈 3개

### 1) Sigma_shift (레벨 변화: “k-시그마 룰”)

-   파라미터: `SIGMA_LEVEL` = k (예: 1.0, 0.5)

**sigma_up:** $$\bar A + k s_A < \bar B - k s_B$$

**sigma_down:** $$\bar A - k s_A > \bar B + k s_B$$

-   **해석**: “원시그마(또는 영쩜오시그마)로 봤을 때 밴드가 아예 벌어졌나?”

### 2) Wilcoxon (분포 중앙/순서 변화)

-   칩 전체 풀링 후 A vs B Wilcoxon rank-sum

-   파라미터: `WILCOX_ALPHA` (강하게, 예: 1e-3)

-   `wilcox_flag` = p \< `WILCOX_ALPHA`

-   **해석**: “평균은 애매한데, 분포가 전체적으로 밀렸나?”

### 3) KS (Radius 기반 공간 성향 변화)

-   Radius(signed) 축 기반으로 공간 특성을 1D로 요약해서 KS로 비교

-   파라미터: `KS_ALPHA` (강하게, 예: 1e-3)

-   `ks_flag` = p \< `KS_ALPHA`

-   **해석**: “좌/우 위치 축에서 패턴(형태)이 바뀌었나?” (센터/엣지 개념이 아니라, 너 실무 Radius 정의에 맞춘 해석)

## results.csv v1.0 최종 스펙 (대소문자 규칙 포함)

### ✅ 컬럼명 규칙

-   소문자 + snake_case 통일
-   msr 값은 원본 이름 그대로(예: `MSR_ABC`, `EDS Hot Bin`, 등) 유지 가능

### ✅ 컬럼 목록(순서 고정)

1.  `msr`
2.  `sigma_up`
3.  `sigma_down`
4.  `wilcox_flag`
5.  `ks_flag`
6.  `mean_ref`
7.  `sd_ref`
8.  `mean_target`
9.  `sd_target`
10. `mean_diff`
11. `median_diff`

### ✅ 컬럼 의미

-   `msr`: MSR 컬럼명(= raw에서 PARTID 이후 열 이름)
-   `sigma_up` / `sigma_down`: Sigma_shift 플래그
-   `wilcox_flag`: Wilcoxon 유의 플래그
-   `ks_flag`: Radius 기반 KS 유의 플래그
-   `mean_a`, `sd_a`, `mean_b`, `sd_b`: 칩 풀링 요약 통계
-   `mean_diff`: `mean_b` - `mean_a`
-   `median_diff`: `median_b` - `median_a`

## 운영 철학(현업용 한 줄)

-   `sigma_shift` = “큰 변화(밴드 분리)”
-   `wilcox` = “조용한 drift”
-   `ks` = “공간 성향 변화(좌/우 축)”

결론은 자동으로 내리지 않고, flag/score로 우선순위를 만들어 사람이 본다.
