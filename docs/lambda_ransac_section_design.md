# Section 3.1–3.2 λ 선정과 RANSAC 연결 설계

## 0. 문서 목적

이 문서는 다음 질문에 대한 실험·논문 서술 원칙을 정리한다.

1. 정규화 파라미터 λ를 전체 visible 점에서 정할지, RANSAC 크기의 소표본에서 정할지
2. 중심 보정(`DisplacementLocal`)을 λ 선정 및 RANSAC과 어떻게 분리할지
3. L-curve, GCV, CD, HD를 어느 단계의 모델에서 계산할지
4. Section 3.1에서 정한 λ를 Section 3.2의 outlier/RANSAC 실험으로 어떻게 넘길지

핵심 결론은 다음과 같다.

> RANSAC에는 소표본 hypothesis와 전체 consensus-inlier refit이라는 서로 다른 추정 문제가 있다. 따라서 Section 3.1에서 `λ_hyp`와 `λ_refit`의 적정 범위를 각각 확인하고, Section 3.2에서는 이를 다시 튜닝하지 않은 채 outlier rejection과 수렴을 검증한다.

---

## 1. λ의 두 역할

현재 solver는 하나의 `lambda`를 provisional hypothesis와 local refit에 모두 전달한다. 그러나 두 단계의 데이터 조건은 본질적으로 다르다.

| 구분 | 역할 | 사용 데이터 | 주요 문제 | 권장 기호 |
|---|---|---|---|---|
| RANSAC hypothesis | 좋은 consensus 후보 생성 | `k × nT` 소표본 | rank 부족, 보간, coverage 편향, 표본 중심 오차 | `λ_hyp` |
| Final refit | 최종 곡면과 중심 추정 | 최종 consensus inlier 전체 | bias–variance, 비관측 영역 extrapolation | `λ_refit` |

현재 회귀는 다음처럼 N-normalized되어 있다.

\[
\min_\beta \frac{1}{N}\|\Phi\beta-\mathbf{1}\|^2
+ \lambda\,\beta^\top W\beta,
\]

즉 구현에서는 normal equation에 `N λ W`를 더하므로 λ의 명목상 크기는 표본 수 N에 덜 민감하다. 하지만 다음 이유로 소표본과 전체점의 최적 λ가 반드시 같지는 않다.

- 소표본의 방향·공간 coverage가 매번 다름
- `k=1.0`이면 28개 계수에 28개 점을 사용하여 OLS가 거의 보간함
- `k<1.0`이면 underdetermined system이 될 수 있음
- 소표본의 설계행렬 condition number가 크게 변함
- 표본 평균이 실제 물체 중심 및 전체 visible centroid와 다를 수 있음
- occlusion이 심할수록 관측 영역이 국소화되어 같은 N에서도 정보량이 감소함

따라서 N-normalization은 λ의 단위를 맞추는 장치이지, 모든 표본 크기와 관측 기하에서 동일한 최적 λ를 보장하는 장치가 아니다.

---

## 2. 전체 Section 흐름

```text
Section 3.1: 무엇을 사용할 것인가
  ├─ 3.1-A  기저 및 조건수 비교
  ├─ 3.1-B  전체점 final-refit λ_refit 선정
  ├─ 3.1-C  소표본 hypothesis λ_hyp 안정성 분석
  ├─ 3.1-D  λ_hyp/λ_refit 단일화 여부 결정
  └─ 3.1-E  선택 λ로 중심 refinement 검증

Section 3.2: 얼마나 반복해야 하며 outlier에 강건한가
  ├─ 3.2-A  λ_hyp/λ_refit 고정
  ├─ 3.2-B  outlier 비율 sweep
  ├─ 3.2-C  omega 갱신 유/무
  ├─ 3.2-D  iteration별 best-model 수렴
  └─ 3.2-E  최적 iteration과 최종 성능 결정
```

Section 3.1은 λ를 정하는 절이고, Section 3.2는 정해진 λ가 RANSAC 안에서 유효한지 검증하는 절이다. Section 3.2에서 outlier 비율마다 GT를 이용해 λ를 다시 최적화하면 두 절의 역할이 무너진다.

---

## 3. Section 3.1-A: 기저와 조건수

### 목적

- Monomial과 Fischer 기저의 수치 조건 비교
- Isotropic ridge와 Sobolev ridge의 구조적 차이 확인

### 원칙

- 같은 clean point cloud 사용
- 같은 차수와 항 수 사용
- 중심과 좌표 스케일 고정
- RANSAC 및 random minimal sampling 사용하지 않음

### 지표

- `cond(Φ'Φ/N)` 또는 eigenvalue spectrum
- effective rank
- coefficient seminorm
- 계산시간

---

## 4. Section 3.1-B: 전체점에서 λ_refit 선정

### 목적

최종 consensus inlier를 모두 사용하는 refit 단계의 λ를 정한다.

### 권장 데이터 흐름

```text
clean partial point cloud
→ 전체 visible 점을 inlier로 사용
→ 중심 좌표 고정
→ 동일한 Φ와 target을 모든 λ에서 공유
→ λ별 최종 ridge 해 계산
→ L-curve/GCV 및 GT 기반 CD/HD 계산
```

### 왜 중심을 고정해야 하는가

고전적인 L-curve는 고정된 `(A,b)`에서 λ만 바뀌는 해의 족을 전제로 한다. λ마다 중심 refinement까지 다시 수행하면 좌표와 설계행렬 Φ가 달라져 순수한 L-curve 해석이 흐려진다.

따라서 λ_refit의 주 L-curve에서는 다음을 고정한다.

- point set
- 중심
- 좌표 스케일
- noise realization
- occlusion/view

### L-curve 축

- y축: 전체 visible 점에서의 residual norm 또는 MSE
- x축: `β'Wβ` 또는 `||Lβ||²`

주의: x축에 λ를 다시 곱한 `λβ'Wβ`를 쓰면 λ=0에서 x=0이 되고, 고전적 solution seminorm의 의미도 바뀐다. 현재 코드처럼 `β'Wβ`를 쓰되 그림 라벨에서 W 또는 L의 정의를 명시한다.

### GT 기반 oracle 지표

- complete clean GT와 양방향 CD
- hidden GT 영역의 recall 방향 거리
- HD95
- 가능하면 F-score와 volumetric IoU

L-curve/GCV는 GT 없이 λ를 고르는 선택 규칙이고, CD/HD는 합성 GT로 그 선택의 타당성을 확인하는 oracle 검증이다. 두 곡선을 같은 것으로 서술하지 않는다.

---

## 5. Section 3.1-C: 소표본에서 λ_hyp 선정

### 목적

RANSAC proposal 크기의 소표본에서 안정적이고 높은 score를 내는 λ를 정한다.

### 실제 solver를 모사해야 하는 부분

```text
전체 clean visible 점 P
→ k × nT개를 추출
→ 추출점 평균으로 전체 좌표를 unbias
→ 추출점에서 provisional β 계산
→ 전체 visible 점에서 residual과 preScore 계산
→ 여러 seed/draw에 대해 반복
```

현재 solver의 provisional 단계와 동일하게 소표본에서 모델을 만들되, 평가는 전체 visible 점에서 해야 한다.

### 소표본 training MSE로 L-curve를 만들면 안 되는 이유

`k=1.0`에서는 28개 계수에 28개 점이므로 OLS training residual이 거의 0이 된다. 이 값은 generalization이나 hidden geometry 품질을 의미하지 않는다. 현재 L-curve에서 `λ=0` MSE가 `10^-20` 이하로 내려간 현상이 이 문제다.

따라서 소표본 단계는 엄밀한 L-curve가 아니라 다음 이름으로 구분한다.

> Proposal stability/generalization vs λ

### 권장 지표

- 전체 visible 점의 `preScore`
- 전체 visible residual 또는 held-out residual
- `P(preScore > threshold)`
- condition number/effective rank
- coefficient seminorm
- CD/HD의 median, IQR
- `P(CD < threshold)`
- spurious component 또는 비정상 등위면 발생률

CD/HD는 소표본 hypothesis의 품질 분포를 진단하는 보조 지표다. 이를 최종 RANSAC 복원 성능으로 보고하지 않는다.

### 샘플링 조건

- 실제 RANSAC과 같은 k 사용
- λ 간에는 동일한 draw를 공유하여 paired 비교
- view/occlusion 간에는 사전에 정의한 동일 view를 사용
- 복원추출 여부를 실제 solver와 맞춤

현재 `lcurve_demo`는 `randperm`으로 비복원추출하지만 실제 solver는 `randsample(..., true, omega)`로 복원추출한다. proposal 실험에서는 이 차이를 의도적으로 통제하거나 실제 solver와 맞춰야 한다.

---

## 6. Section 3.1-D: λ를 하나로 쓸지 두 개로 쓸지

### 1단계: 두 basin 확인

- 전체점 oracle/L-curve에서 안정적인 `λ_refit` 범위
- 소표본 score/CD 분포에서 안정적인 `λ_hyp` 범위

### 두 범위가 겹치는 경우

하나의 λ를 사용한다.

- 두 실험의 안정적인 plateau 교집합에서 선택
- 단일 조건의 날카로운 최소값보다 여러 view/occlusion에서 안정적인 중앙 영역 선호
- 장점: 파라미터 수가 적고 논문·배치 설명이 단순함

### 두 범위가 겹치지 않는 경우

solver 파라미터를 분리한다.

```matlab
ransacPar.lambdaHyp
ransacPar.lambdaRefit
```

- provisional `betaTmp`: `lambdaHyp`
- consensus inlier의 `DisplacementLocal` 및 최종 β: `lambdaRefit`

소표본 hypothesis와 최종 inlier refit은 데이터 조건과 목적이 다르므로 두 λ를 사용하는 것은 방법론적으로 정당하다. 단, 데이터가 두 λ 분리 필요성을 명확히 보여줄 때만 도입한다.

---

## 7. Section 3.1-E: 중심 refinement 검증

### RANSAC과 분리할 부분

현재 `DisplacementLocal`은 RANSAC solver 내부 nested function이지만, 중심 보정 자체는 outlier rejection과 다른 기능이다. 다음과 같은 독립 함수로 분리하는 것이 바람직하다.

```matlab
[center, beta, history] = refineCenterAndFit(points, Funcs, Grads, lambda, opts)
```

권장 출력:

- 최종 중심 및 β
- iteration별 center shift
- residual history
- convergence flag
- 사용 iteration 수

권장 종료 조건:

- `maxIters`
- `tolCenter`
- 필요하면 `tolObjective`

현재처럼 항상 4회만 수행하는 방식은 알고리즘 정의에는 간단하지만, 실제 수렴 여부를 보여주지 못한다.

### 중심 보정 평가에서 전체점을 써야 하는 이유

중심 refinement 자체를 평가할 때는 전체 clean visible inlier를 사용한다. `k×nT`개만 사용하면 다음이 섞인다.

- occlusion으로 인한 visible-centroid bias
- random subsampling으로 인한 centroid noise
- 중심 refinement 알고리즘 오차

소표본 중심화는 Section 3.1-C의 RANSAC hypothesis 안정성에서 평가하고, 최종 중심 refinement는 전체 inlier에서 평가한다.

### 의미 있는 중심 강건성 조건

모든 점에 동일한 global translation을 더하는 실험은 현재 mean-centering 구조상 거의 자명한 translation equivariance 확인에 그친다. 더 중요한 조건은 다음이다.

- occlusion으로 visible centroid가 GT center에서 벗어나는 정도
- view 방향
- 비균일한 표면 밀도
- noise 수준
- 부분적으로 한쪽 면만 관측되는 경우

권장 지표:

- 중심 오차 `||ĉ-c_GT||`
- refinement 전후 CD/HD
- 수렴 iteration
- 실패율

---

## 8. Section 3.2: Outlier와 RANSAC 수렴

### 기본 원칙

Section 3.2에서는 Section 3.1에서 선택한 `λ_hyp`, `λ_refit`을 고정한다. 목적은 λ를 다시 고르는 것이 아니라 다음을 평가하는 것이다.

- outlier rejection
- adaptive omega의 효과
- iteration에 따른 best-model 품질
- 충분한 고정 반복 횟수
- clean에서 고른 λ의 outlier 환경 이식성

### 권장 RANSAC 흐름

```text
k × nT 소표본
→ λ_hyp로 provisional hypothesis
→ 전체 점에서 score/inlier mask
→ local optimization 후보
→ best consensus 선택
→ best inlier 전체
→ λ_refit으로 중심 refinement + 최종 β
→ 최종 CD/HD/F1
```

최종 모델을 minimal sample만으로 유지하지 않는다. 소표본은 consensus를 찾기 위한 proposal이고, 보고할 최종 표면은 전체 consensus inlier로 refit한 결과다.

### λ 비교 범위

주 결과:

- `λ_hyp`, `λ_refit`: Section 3.1 결정값

ablation:

- OLS: `(0,0)`
- 선택 λ
- 선택 λ의 10배 또는 1/10배

outlier 비율마다 GT로 새 최적 λ를 선택하지 않는다. 만약 선택 λ가 특정 outlier 조건에서 실패한다면, per-condition oracle retuning보다 post-consensus L-curve 같은 고정된 선택 규칙을 검토한다.

### CD/HD 계산 대상

- 모든 폐기 hypothesis에 CD/HD를 계산하여 평균하지 않음
- iteration별로는 그 시점의 running-best model에서 계산
- 최종 결과는 best consensus inlier refit 모델에서 계산
- F1/precision/recall은 inlier classification 성능
- CD/HD는 최종 geometry 성능

### 반복 종료

현재 MC 코드는 `freezeIter=true`이므로 confidence 기반 조기 종료 없이 `maxIter=2000`을 고정 실행한다. 이는 λ 효과와 미수렴 효과를 분리하는 통제 실험에는 적합하다. 다만 다음 saturation 확인이 필요하다.

> 2000회와 4000회에서 최종 CD/F1과 λ 순위가 실질적으로 변하지 않는가?

현재 local refinement도 `locIters=4` 고정이며 tolerance 기반 수렴은 아니다. RANSAC 반복 수렴과 local center refinement 수렴을 구분해 보고한다.

---

## 9. L-curve, GCV, CD/HD의 단계별 역할

| 단계 | L-curve | GCV | CD/HD | 핵심 출력 |
|---|---|---|---|---|
| 전체점 final refit | 주 선택 도구 | 보조 선택 도구 | oracle 검증 | `λ_refit` |
| 소표본 hypothesis | 사용하지 않거나 참고만 | held-out 방식 권장 | 품질 분포 진단 | `λ_hyp` |
| 중심 refinement | λ 선택 후 고정 | 불필요 | 전후 성능 비교 | 중심 강건성 |
| RANSAC iteration | 사용하지 않음 | 사용하지 않음 | running-best에서 계산 | 수렴 곡선 |
| 최종 RANSAC | 필요 시 final inlier 진단 | 선택적 | 주 최종 지표 | end-to-end 성능 |

### CD/HD 해석 주의

- 합성 실험의 GT는 noisy/occluded 입력이 아니라 complete clean surface를 사용
- hidden-region metric과 complete-surface metric을 구분
- 서로 다른 occlusion 비율은 hidden target set 자체가 달라지므로 absolute `CD_occ` 크기의 단순 비교에 주의
- occlusion 간 결론은 같은 view의 `ridge − OLS` 또는 상대 개선량을 우선 사용
- HD95는 spurious lobe 억제를 잘 반영하지만 CD 최적 λ보다 강한 λ를 선호할 수 있음

---

## 10. 통계 및 재현성 원칙

1. 동일 `(shape, view, occlusion, draw)`에서 λ를 paired 비교한다.
2. λ가 바뀌어도 noise, view, 소표본 index를 공유한다.
3. occlusion별 비교에는 동일 view와 동일 noisy base cloud를 사용한다.
4. 소표본 반복은 독립 물체 수로 세지 않고 해당 view 내부 반복으로 취급한다.
5. 평균±표준편차와 함께 median/IQR 또는 성공률을 보고한다.
6. corner는 평균 곡선 하나뿐 아니라 view별 corner λ 분포도 확인한다.
7. 최종 λ를 같은 test view에서 선택하고 평가하지 않도록 조건 또는 shape를 train/validation으로 분리하는 것이 이상적이다.

---

## 11. 현재 코드에서 바로 유의할 점

### `lcurve_demo.m`

- 실제 RANSAC을 호출하지 않음
- `k×28`개 random subset에 직접 ridge를 수행함
- 소표본 자체의 MSE를 사용하여 `k=1.0, λ=0`에서 보간 잔차가 퇴화함
- MSE/penalty는 draw 10회의 중앙값이지만 CD/HD는 첫 draw 하나만 사용함
- 현재 결과는 final-model L-curve라기보다 proposal-level regularization 진단에 가까움

### `PoliNavigationSolver3_3_FischerRansac_MC.m`

- provisional fit은 복원추출 `randsample(..., true, omega)` 사용
- λ가 provisional과 local refit에 동일하게 들어감
- 최종 best model은 local optimization을 수행한 후보 중 최고 postScore 모델
- 별도의 마지막 full-inlier refit을 명시적으로 한 번 더 수행하지 않음
- `freezeIter=true`이면 고정 반복이며 convergence stop이 아님
- local refinement는 tolerance 없이 고정 `locIters`

### `test_ransac_ridge_montecarlo_image.m`

- 현재 활성 λ 목록에는 정확한 `λ=0` 대신 `1e-6`이 들어 있음
- OLS baseline을 주장하려면 실제 `λ=0`을 포함해야 함
- 주석은 “모든 적응 메커니즘 동결”이라고 되어 있으나 현재 설정은 `freezeW=false`, `freezeOmega=false`임
- k struct 초기값은 loop에서 덮어쓰므로 실제 값은 loop의 `k_list`임

---

## 12. 실행 우선순위

### 1단계: 함수 분리

- `fitSobolevRidge(points, lambda, ...)`
- `refineCenterAndFit(points, lambda, opts)`
- RANSAC solver는 두 함수를 호출하도록 정리
- 필요하면 `lambdaHyp`와 `lambdaRefit` 분리

### 2단계: 전체점 λ_refit 실험

- 중심 고정
- full visible inlier
- L-curve/GCV/CD/HD
- view/occlusion별 basin 확인

### 3단계: 소표본 λ_hyp 실험

- 실제 RANSAC sampling과 동일 조건
- 전체 visible score로 평가
- k/view/occlusion별 성공률 및 CD/HD 분포

### 4단계: 단일 λ 여부 결정

- 두 basin이 겹치면 단일 λ
- 겹치지 않으면 `lambdaHyp`, `lambdaRefit` 분리

### 5단계: 중심 refinement

- 전체 visible inlier
- 선택 `λ_refit`
- occlusion/view별 중심 오차와 수렴 분석

### 6단계: Section 3.2 RANSAC

- 선택 λ 고정
- outlier 및 iteration sweep
- omega 갱신 ablation
- running-best 및 final-refit CD/HD/F1
- 2× iteration saturation 확인

---

## 13. 논문용 핵심 서술

> Regularization was characterized separately for the hypothesis-generation and consensus-refinement stages of RANSAC. The hypothesis penalty was selected from repeated minimal-sample fits evaluated on the full clean observation, whereas the refinement penalty was selected using an L-curve and GCV on a fixed full-inlier system and validated against geometric oracle metrics. These parameters were then fixed during the outlier and convergence experiments, preventing regularization tuning from being confounded with RANSAC robustness.

단일 λ를 사용할 수 있는 경우에는 다음처럼 축약한다.

> The normalized Sobolev penalty exhibited an overlapping stable basin in both minimal-sample hypothesis generation and full-inlier refinement. A single λ was therefore fixed for all subsequent RANSAC and outlier experiments.

