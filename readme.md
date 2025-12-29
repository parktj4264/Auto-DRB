# DRB rev1 (Final Version)

## 1. 개요 및 목표

**A공정(REF) → B공정(TARGET) 변경 시, 유의미한 변동이 있는 MSR을 정량적으로 평가하고 랭킹화한다.**

-   본 알고리즘은 **판정(O/X)**을 내리는 것이 아니라, 엔지니어가 불량의 심각도(Magnitude)를 인지할 수 있도록 **점수(Score)와 방향(Direction)**을 제공하는 "리스트업 필터(Screening Filter)"다.
-   기존의 P-value 기반 검정을 지양하고, **Effect Size(효과 크기)**를 통해 실질적인 차이를 측정한다.

------------------------------------------------------------------------

## 2. 데이터 전제 (Input Spec)

### `raw.csv` (Chip Level Data)

-   **Key**: `ROOTID`
-   **필수 컬럼**: `PARTID`, `Radius`
-   **MSR 컬럼**: `PARTID` 컬럼 바로 다음 열부터 끝까지 전부
-   **📌 Radius 주의사항**:
    -   실무의 `Radius`는 이름만 Radius고, 실제 의미는 **Signed Position (좌/우 위치)**로 가정한다.
    -   범위: 대략 -150 \~ 150 (실수형)

### `ROOTID.csv` (Wafer Level Meta)

-   **필수 컬럼**: `ROOTID`, `GROUP`
-   **GROUP 값**: `REF` (기준), `TARGET` (비교군)

------------------------------------------------------------------------

## 3. 분석 렌즈 (3 Scores + Direction)

DRB rev1은 3가지 관점의 **점수(Score)**와 1가지 **방향성(Flag)**을 산출한다.

### 1) Level Shift: `sigma_score` (Glass's Delta)

Reference 공정의 산포를 기준으로 Target 공정의 평균이 얼마나 이동했는지 측정한다. \* **파라미터**: 없음 (수식 자체로 정의됨) \* **수식**: $$\text{score} = \frac{\mu_{target} - \mu_{ref}}{\sigma_{ref}}$$ \* **해석**: "Ref 산포 대비 몇 시그마만큼 튀었는가?" (+3.0이면 3시그마 상향 돌파)

### 2) Global Drift: `cliffs_delta` (Stochastic Dominance)

전체 분포의 서열 관계를 통해 확률적 우세(쏠림)를 측정한다. \* **구현**: Wilcoxon Rank Sum Test의 $W$ 통계량 역산 \* **수식**: $$\delta = \frac{2W - n_1 n_2}{n_1 n_2} \quad (\text{Target} > \text{Ref} \text{ 일 때 양수})$$ \* **범위**: -1.0 \~ 1.0 \* **해석**: "랜덤하게 칩을 뽑았을 때 Target이 클 확률이 얼마나 압도적인가?" (절댓값이 클수록 분포 분리 심각)

### 3) Spatial Shape: `ws_spatial` (Profile $L_1$ Deviation)

Radius 축(Shot map trend)에 따른 공간적 프로파일의 **물리적 괴리도(Physical Discrepancy)**를 측정한다.

-   **구현 (Implementation)**

    1.  Radius를 $K$개의 구간(Bin)으로 분할하여 국소 평균(Local Mean) 벡터를 생성한다.
    2.  Ref와 Target 프로파일 간의 $L_1$ 거리 (Mean Absolute Difference)를 계산한다. \> *Note: 반도체 공정 특성상 칩의 수평 이동(Transport)보다는 **위치별 불량 심도(Vertical Gap)**가 중요하므로, 연산 효율과 공정 적합성을 고려하여 Binning* $L_1$ 방식으로 형상 차이를 근사함.

-   **수식** $$
    \text{score} = \frac{1}{K} \sum_{i=1}^{K} \left| \text{Profile}_{ref}(i) - \text{Profile}_{target}(i) \right|
    $$

-   **해석**: **"공간적 트렌드(Shape)가 물리적으로 얼마나 망가졌는가?"**

    -   값이 **0에 가까울수록**: 두 공정의 Shot Map 형태가 완벽하게 일치함.
    -   값이 **클수록**: Center/Edge 경향성이나 국소적인 튐(Spike) 현상이 발생하여 모양이 달라짐.

### 4) Direction Flag: `direction` (Up/Down/Stable)

엔지니어의 직관적 판단을 돕기 위해 변화의 방향을 명시한다. \* **파라미터**: `DIR_SIGMA_THRESHOLD` (예: 1.0, 0.5 등 사용자가 민감도 조절 가능) \* **로직**: \* `sigma_score` $\ge$ `DIR_SIGMA_THRESHOLD` $\rightarrow$ **"Up"** \* `sigma_score` $\le$ -`DIR_SIGMA_THRESHOLD` $\rightarrow$ **"Down"** \* 그 외 $\rightarrow$ **"Stable"**

------------------------------------------------------------------------

## 4. 결과물 스펙 (`results.csv`)

### ✅ 컬럼명 규칙

-   소문자 + snake_case 통일
-   msr 값은 원본 이름 그대로 유지

### ✅ 컬럼 목록 (순서 고정)

| 순서 | 컬럼명 | 설명 | 비고 |
|:-----------------|:-----------------|:-----------------|:-----------------|
| 1 | **msr** | MSR 컬럼명 | Key |
| 2 | **direction** | 변화 방향 (Up/Down/Stable) | **직관적 필터** |
| 3 | **sigma_score** | Glass's Delta (레벨 변동 점수) | **주요 Sorting 기준** |
| 4 | **cliffs_delta** | Cliff's Delta (분포 쏠림 점수) | 보조 지표 (Shift) |
| 5 | **ws_spatial** | Spatial Wasserstein (모양 변화 점수) | 보조 지표 (Shape) |
| 6 | **mean_ref** | Reference 평균 |  |
| 7 | **sd_ref** | Reference 표준편차 |  |
| 8 | **mean_target** | Target 평균 |  |
| 9 | **sd_target** | Target 표준편차 |  |

------------------------------------------------------------------------

## 5. 운영 철학 (Engineering Guide)

결과는 **`sigma_score` 절대값 내림차순**으로 정렬하여 보는 것을 권장한다.

| 지표             | 현업용 해석 (한 줄 요약)                                |
|:-----------------------------------|:-----------------------------------|
| **sigma_score**  | "기존 공정(Ref) 관리폭 대비 **몇 배(Sigma)나 튀었나?**" |
| **cliffs_delta** | "분포가 얼마나 **한쪽으로 쏠렸나?**" (-1\~1)            |
| **ws_spatial**   | "맵(Map)의 **생김새가 얼마나 다른가?**"                 |
| **direction**    | "그래서 **열화(Up/Down)야 개선이야?**"                  |
