# 실험 계획 및 논문 구성 정리

**대상 논문:** Even-Degree Homogeneous Polynomial Modeling for 3D Obstacle Representation (HERO)
**투고 타겟:** T-AES Navigation 영역 / Advances in Space Research
**작성일:** 2026년 5월

---

## 0. 핵심 포지셔닝 (한 줄)

"3D LiDAR 점군 → Fischer 분해 기반 짝수차 동차다항식 음함수 피팅"을 비협조 우주 도메인에서 최초로 완결하는 피팅 방법론 논문.  
제어 통합(CBF/MPC)은 2편에서 다루며, 본 논문은 **C1~C4 기여 + ISS 정적 멀티뷰 피팅** 범위로 한정한다.

**C4 포지셔닝 주의:** 그래프 Laplacian 기반 스펙트럴 클러스터링 자체는 기존 기법이므로 독립 기여로 내세우지 않는다.  
"HERO 파이프라인의 전처리 단계로 채택하여 다물체 우주 구조물로 자연스럽게 확장됨을 보인다"는 방식으로 포지셔닝.

---

## 1. Section 3 실험 구조

### 논리 흐름

```
3.1 무엇을 쓸 것인가   → 기저 선택 + λ 선정
3.2 얼마나 반복할 것인가 → 수렴 분석 + 최적 iter
3.3 실제 환경에서 되는가 → ISS 멀티뷰 시뮬레이션
      ├── 3.3.1 단일 물체 (CubeSat 등)
      └── 3.3.2 다물체 확장 (ISS) — 스펙트럴 클러스터링 전처리 포함
```

---

### 3.1 Basis Selection and Regularization

**목적:** Fischer 기저의 조건수 이점과 Sobolev Ridge의 정당성을 정량적으로 확립.  
이 절 끝에서 "이후 모든 실험은 Fischer + Sobolev, λ = 최적값"을 한 줄로 고정.

#### 비교 구도 (공정 비교 설계)

| 비교군 | 기저 | 차수 | 항 수 | Ridge 방식 |
|---|---|---|---|---|
| Baseline | Monomial (동차, homogeneous) | 6 | 28 | Isotropic (λI, 전 항 동일) |
| **HERO** | Fischer (동차 구면조화) | 6 | 28 | Sobolev H¹ (level별 차등) |

**동차 단항식 사용 이유:** 비동차 기저는 항 수와 표현력이 달라지므로 "기저 선택" 효과를 분리 불가. 동차로 고정해야 항 수가 동일하여 공정 비교.

**Sobolev vs Isotropic 이유:** 단항식 기저는 자연스러운 level 구조가 없으므로 차등 페널티가 의미 없음. 각 기저에 가장 적합한 페널티를 각각 부여하는 것이 공정한 비교.

#### (a) 조건수 비교

- **실험:** 동일 점군(N=10,000, HexagonPrism, eps=0%)에 두 기저 적용
- **지표:** κ(Φ'Φ/N) = max(eig)/min(eig), eigenvalue 분포
- **Figure:** eigenvalue 분포 세미로그 plot (두 기저 겹쳐 그리기)
- **기대 결과:** Fischer 기저의 κ가 단항식 대비 수 배~수십 배 작음

#### (b) λ 선정 — L-curve

**핵심 아이디어:** MSE는 λ가 클수록 단조 증가 → MSE 단독으로는 최적 λ 결정 불가.  
L-curve (MSE vs ‖β‖²)의 elbow가 bias-variance 최적 지점.

```
MSE ↑
 │·  (λ→0: β 발산, MSE 낮음)
 │ ·
 │  ·
 │   · ← elbow = 최적 λ
 │     ·····  (λ 큼: MSE 급등)
 └──────────── ‖β‖² →
```

- **λ sweep 범위:** [0, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.3, 1]
- **Figure 1:** L-curve (두 Ridge 방식 × 두 기저 → 4개 곡선)
- **Figure 2:** λ vs Chamfer(occluded) 곡선 (L-curve elbow와 일치함을 검증)
- **GCV 보조 검증:** Generalized Cross-Validation score = MSE / (1 - tr(H)/N)², hat matrix H는 ridge에서 closed-form

#### (c) Sobolev vs Isotropic Ridge 비교

동일 λ 예산에서, 두 Ridge 방식의 표면 품질 차이를 보임.

| 지표 | Monomial + Iso | Fischer + Iso | **Fischer + Sobolev** |
|---|---|---|---|
| κ(Φ'Φ/N) | 큼 | **작음** | **작음** |
| ‖β‖² | 큼 | 중간 | **작음** |
| Chamfer(visible) | 기준 | ↓ | ↓ |
| **Chamfer(occluded)** | 기준 | ↓↓ | **↓↓↓** |
| 계산시간(ms) | 기준 | ≈ 같음 | ≈ 같음 |

**Chamfer(occluded)가 이 논문의 킬러 지표:** 동일 관측 데이터로 피팅했을 때 보이지 않는 영역의 표면 품질 → 우주 부분 관측 상황에서의 실용성.

---

### 3.2 Convergence Analysis and Optimal Iteration

**목적:** RANSAC의 ω(오메가) 갱신 메커니즘 효과와 최적 iter 수를 결정.

**고정 조건:** Fischer + Sobolev, λ = 3.1절에서 결정된 최적값.

#### (a) ω 갱신 유/무 비교

`ransacPar.sim.freezeOmega` 플래그로 제어 가능 (코드에 이미 있음).

- **Case A (freezeOmega=true):** 균일 샘플링, omega 갱신 없음
- **Case B (freezeOmega=false):** Gaussian 소프트 스코어 기반 점진적 재가중

```
ω 갱신 있는 경우: 이상치 점의 샘플링 확률이 반복될수록 감소
                  → 후기 iter에서 가설 품질 향상
```

- **Figure:** iter vs F1 (Case A vs B 비교, eps=20% 조건)
- **지표:** 동일 F1 달성에 필요한 iter 수 비교 (ω 갱신이 수렴 속도를 몇 배 개선하는가)

#### (b) iter vs 피팅 품질 (λ 종류별)

여러 λ에 대해 iter vs Chamfer/F1을 겹쳐 그림.  
**메시지:** 최적 λ에서 plateau가 가장 낮고 수렴이 빠름.

```
Chamfer ↓
   │\  (λ 너무 작음: plateau 높음)
   │ \___(λ 최적: plateau 낮음, 빠른 수렴)
   │     \_____(λ 너무 큼: plateau 낮지만 수렴 느림)
   └──────────── iter →
```

- **λ 3종:** 0 (OLS), 최적값, 최적값×10
- **지표:** Chamfer(occ), F1, postScore vs iter
- **Figure:** 3×1 subplot (Chamfer / F1 / postScore), λ별 색 구분

#### (c) 최적 iter 결론

- **기준:** Chamfer가 plateau의 95%에 도달하는 iter = 최적 iter
- **Figure:** Chamfer vs iter + 계산시간 누적 → "X회 iter = Y ms, 품질 Z Chamfer"
- **표:** eps별 최적 iter 요약 (0% / 10% / 20% / 40%)

---

### 3.3 ISS Multi-view Simulation

**목적:** 실제 우주 도메인 데이터에서 HERO가 동작함을 보임.  
**범위 제한 (Paper 1):** 정적 멀티뷰 피팅. CW 동역학 궤적 시뮬레이션은 Paper 2.

---

#### 3.3.1 단일 물체 — CubeSat / 단순 위성 형상

ISS로 넘어가기 전, 단일 연결 형상에서 멀티뷰 피팅이 동작함을 보이는 준비 단계.

**구현 방식:**

```
CubeSat / 단순 위성 CAD 모델
    │
    ├── 뷰포인트 k개 (k = 1, 2, 3, 4, 6)
    │     └── 레이캐스팅 → 부분 점군 + 노이즈 + 이상치(eps=20%)
    │
    └── 순차 누적 피팅
          └── Chamfer(occ, GT) 추적
```

**보여줄 것:**
- Figure: k(뷰 수) vs Chamfer(occ) 감소 곡선
- Figure: 등위면 스냅샷 3장 (k=1, 3, 6)
- Table: k별 Chamfer(vis), Chamfer(occ), 계산시간

---

#### 3.3.2 다물체 확장 — ISS (스펙트럴 클러스터링 전처리)

**왜 필요한가:**  
ISS는 모듈·솔라패널·트러스가 위상학적으로 분리된 구조. 단일 다항식으로는 표현 불가.  
스펙트럴 클러스터링으로 컴포넌트를 분리한 후, 각 클러스터에 HERO를 독립 적용.

**포지셔닝:** 클러스터링 자체가 기여가 아니라, "HERO 파이프라인이 전처리와 결합하여 복잡 다물체 구조에 적용 가능함"이 기여.

**구현 흐름:**

```
ISS 3D 모델 (공개 CAD / NASA dataset)
    │
    ├── 멀티뷰 레이캐스팅 → 누적 점군 P
    │
    ├── [전처리] 스펙트럴 클러스터링
    │     ├── KNN 그래프 구성 (k_nn = 10~20)
    │     │     └── 거리 + 법선 유사도 가중치
    │     ├── 정규화 Laplacian L = D^{-1/2} W D^{-1/2}
    │     ├── 하위 고유벡터 (nc개) → k-means
    │     └── 출력: 클러스터 레이블 {C_1, ..., C_m}
    │
    └── 각 클러스터 C_i에 독립적으로 HERO 피팅
          └── {β_1, ..., β_m} 출력
```

**비교 실험 (Cutting의 효과를 보임):**

| 방법 | Chamfer(occ) | 설명 |
|---|---|---|
| 단일 다항식 (cutting 없음) | 높음 | ISS 전체를 하나의 f=1로 근사 → 부정확 |
| **HERO + 스펙트럴 cutting** | **낮음** | 컴포넌트별 별도 피팅 → 정밀 표현 |

**보여줄 것:**
- Figure 1: ISS 점군 + 클러스터 분리 결과 (컬러 레이블 시각화)
- Figure 2: 컴포넌트별 등위면 시각화 (각각 다른 색, 합성 뷰)
- Figure 3: Chamfer(occ) 비교 — 단일 다항식 vs HERO+cutting
- Table: 클러스터 수, 각 클러스터 Chamfer(occ), 전체 Chamfer(occ), 계산시간

**클러스터 수 결정:**  
정규화 Laplacian 고유값 갭(spectral gap)으로 자동 결정하거나,  
ISS 주요 구성 요소 수 (모듈 본체 1~2, 솔라패널 어레이 좌/우, 트러스 등)에 맞춰 수동 설정.

**주의사항:**
- CW 동역학 언급은 하지 않음 (Conclusion에서 "향후 연구" 한 줄)
- 클러스터링 파라미터(k_nn, nc) 민감도 분석을 부록 또는 간단한 ablation으로 포함 검토

---

## 2. 구현 방식

### 2.1 비교군 설정

```matlab
% Fischer 기저 (HERO)
TermsF = homogeneFischerTerms(order);          % order=6 → 28항
[FuncsF, GradsF] = makeFuncsGradsStack(TermsF);

% Monomial 동차 기저 (비교군)
% homogeneTerm.m 활용 또는 직접 생성
% x^a * y^b * z^c, a+b+c = order
TermsM = homogeneTerm(order);                  % 동차 단항식
[FuncsM, GradsM] = makeFuncsGradsStack(TermsM);

% Ridge 방식
% Sobolev: makeRidgeDiag(order, lambda, s)
% Isotropic: lambda * eye(nT)
```

### 2.2 Occlusion 생성 방법

```matlab
% 뷰포인트 방향 벡터 v = [vx, vy, vz] (단위벡터)
% 해당 방향에서 보이는 점: 법선 n과 뷰 방향이 반대 방향인 점
% (dot(n_i, v) < 0 인 점만 선택)

function pts_visible = applyOcclusion(pts_full, normals, view_dir)
    dots = normals * view_dir(:);
    pts_visible = pts_full(dots < 0, :);
end

% GT(전체) vs Visible(관측) 분리
% Chamfer(visible): pts_visible ↔ 등위면 샘플
% Chamfer(occluded): pts_full(dots >= 0) ↔ 등위면 샘플
```

### 2.3 Chamfer Distance 계산

```matlab
function cd = chamferDistance(pts_A, pts_B)
    % pts_A → pts_B: 각 점에서 최근접 거리의 평균
    % pts_B → pts_A: 반대 방향
    % Chamfer = 0.5*(mean(d_AB) + mean(d_BA))
    
    d_AB = mean(min(pdist2(pts_A, pts_B), [], 2));
    d_BA = mean(min(pdist2(pts_B, pts_A), [], 2));
    cd = 0.5*(d_AB + d_BA);
end

% 등위면 샘플 pts_B 생성: fimplicit3 또는 isosurface 후 샘플링
```

### 2.4 조건수 계산

```matlab
% 이미 test_ransac_montecarlo_eordnT_iter_V3.m에 있음
X_F = calculateFourthOrder(PPm_use, FuncsF);  % Fischer
X_M = calculateFourthOrder(PPm_use, FuncsM);  % Monomial

ev_F = sort(eig(X_F'*X_F / N), 'descend');
ev_M = sort(eig(X_M'*X_M / N), 'descend');

kappa_F = ev_F(1)/ev_F(end);
kappa_M = ev_M(1)/ev_M(end);
```

### 2.5 스펙트럴 클러스터링 (다물체 전처리)

```matlab
function labels = spectralCutting(pts, k_nn, nc)
% pts  : N×3 점군
% k_nn : KNN 그래프 이웃 수 (권장 10~20)
% nc   : 클러스터 수 (ISS 구성요소 수 또는 spectral gap으로 결정)

    N = size(pts, 1);

    % 1) 법선 추정 (pcnormals 또는 직접 PCA)
    pc = pointCloud(pts);
    normals = pcnormals(pc, k_nn);   % N×3

    % 2) KNN 그래프 가중치 행렬
    [idx, dist] = knnsearch(pts, pts, 'K', k_nn+1);
    idx  = idx(:, 2:end);    % 자기 자신 제외
    dist = dist(:, 2:end);

    sigma_d = median(dist(:));
    sigma_n = 0.5;           % 법선 유사도 가중치 스케일

    W = sparse(N, N);
    for i = 1:N
        for j_idx = 1:k_nn
            j = idx(i, j_idx);
            w_d = exp(-dist(i,j_idx)^2 / (2*sigma_d^2));
            w_n = exp(-(1 - abs(dot(normals(i,:), normals(j,:))))^2 / (2*sigma_n^2));
            w   = w_d * w_n;
            W(i,j) = w;
            W(j,i) = w;
        end
    end

    % 3) 정규화 Laplacian
    D_vec = sum(W, 2);
    D_inv_sqrt = spdiags(1./sqrt(max(D_vec, 1e-12)), 0, N, N);
    L_sym = speye(N) - D_inv_sqrt * W * D_inv_sqrt;

    % 4) 하위 nc개 고유벡터
    [V, ~] = eigs(L_sym, nc, 'smallestabs');  % N×nc

    % 5) 행 정규화 후 k-means
    V_norm = V ./ max(sqrt(sum(V.^2, 2)), 1e-12);
    labels = kmeans(V_norm, nc, 'Replicates', 10, 'MaxIter', 300);
end

% 사용 예시 (ISS 점군)
nc     = 4;   % ISS 주요 구성요소 수 (spectral gap 보고 결정)
labels = spectralCutting(PPm_iss, 15, nc);

for c = 1:nc
    pts_c = PPm_iss(labels == c, :);
    [~, Beta_c, ~, ~] = PoliNavigationSolver3_3_FischerRansac_MC(...
        0, pts_c, order, nT, Funcs1, Grads1, ransacPar);
    % Beta_c: c번째 컴포넌트의 HERO 계수
end
```

**spectral gap으로 nc 자동 결정 (선택):**

```matlab
[~, eigvals] = eigs(L_sym, 10, 'smallestabs');
ev = sort(diag(eigvals));
gaps = diff(ev);
[~, nc_auto] = max(gaps);   % gap이 가장 큰 위치 = 클러스터 수
```

### 2.6 L-curve 생성

```matlab
lambda_list = [0, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.3, 1];
mse_list = zeros(size(lambda_list));
betanorm_list = zeros(size(lambda_list));

for i = 1:numel(lambda_list)
    lam = lambda_list(i);
    beta = regressionFourthOrder(pts_fit, FuncsF, lam, order, 1);
    r = abs(FuncsF(pts_all) * beta - 1);
    mse_list(i) = mean(r.^2);
    betanorm_list(i) = beta'*beta;
end

% L-curve plot
figure; loglog(betanorm_list, mse_list, 'o-');
xlabel('||β||²'); ylabel('MSE');
title('L-curve');
```

---

## 3. 지표 전체 목록

### 3.1 표면 복원 품질

| 지표 | 기호 | 설명 | 절 |
|---|---|---|---|
| **Chamfer Distance (visible)** | CD_vis | 관측 점군 ↔ 등위면 샘플, 양방향 평균 | 3.1, 3.2, 3.3 |
| **Chamfer Distance (occluded)** | CD_occ | 비관측 GT 점군 ↔ 등위면 샘플 | **3.1, 3.3 (킬러 지표)** |
| Hausdorff Distance | HD | 양방향 최근접 최댓값 (최악 오차) | 3.3 보조 |

### 3.2 수치 안정성 / 정규화 효과

| 지표 | 기호 | 설명 | 절 |
|---|---|---|---|
| **조건수** | κ | max(eig)/min(eig) of Φ'Φ/N | 3.1 |
| **‖β‖²** | betaNorm | 계수 L2 노름 제곱, ridge 억제 정도 | 3.1, 3.2 |
| MSE (algebraic) | MSE | mean((Φβ-1)²), ridge 대가 정량화 | 3.1 (L-curve 축) |

MSE 포지셔닝 주의: **"낮을수록 좋은 지표"가 아니라 "ridge 대가의 척도"**. CD_occ가 개선되는 trade-off를 보여주는 용도.

### 3.3 이상치 강건성 (RANSAC)

| 지표 | 기호 | 설명 | 절 |
|---|---|---|---|
| **F1** | F1 | 2·P·R/(P+R), 인라이어 분류 종합 | 3.2 |
| Precision | P | TP/(TP+FP) | 3.2 |
| Recall | R | TP/(TP+FN) | 3.2 |

### 3.4 수렴 / 효율

| 지표 | 기호 | 설명 | 절 |
|---|---|---|---|
| preScore (Sj) | - | raw 가설 Gaussian 소프트 스코어 | 3.2 |
| postScore (SjLoc) | - | local 최적화 후 스코어 | 3.2 |
| inlierR | - | 인라이어 비율 | 3.2 |
| 계산시간 | t [ms] | 피팅 총 소요시간 | 3.2, 3.3 |

### 3.5 ISS 시나리오 전용

| 지표 | 기호 | 설명 | 절 |
|---|---|---|---|
| CD_occ(k) | - | 누적 뷰 k에서의 occluded Chamfer | 3.3.1 |
| 커버리지(k) | cov | k개 뷰에서 관측된 표면 비율 (%) | 3.3.1 |

### 3.6 다물체 분리 전용 (3.3.2)

| 지표 | 기호 | 설명 |
|---|---|---|
| **CD_occ (단일 다항식)** | - | Cutting 없이 ISS 전체 하나로 피팅한 결과 |
| **CD_occ (HERO + cutting)** | - | 컴포넌트별 분리 후 피팅 결과 — 차이가 메시지 |
| Precision / Recall (클러스터) | - | GT 컴포넌트 레이블과 비교한 클러스터링 정확도 |
| ARI (Adjusted Rand Index) | - | 클러스터링 품질 종합 지표 |
| 컴포넌트별 CD_occ | - | 각 클러스터의 개별 피팅 품질 |
| 계산시간 (클러스터링 + 피팅) | - | 전처리 오버헤드 정량화 |

---

## 4. 지표 선정 근거 요약

Ridge를 추가할수록 MSE는 증가 (bias 증가)하지만 ‖β‖²는 감소 (variance 감소).  
이 trade-off를 시각화하는 것이 L-curve이고, 실제로 유익한 경우는 **보이지 않는 영역에서 표면 품질이 개선될 때**이다.

따라서:
- **λ 선정 기준:** L-curve elbow + CD_occ 최소화 지점의 일치 확인
- **기저 선정 기준:** 동일 λ 예산 아래 κ와 CD_occ 비교
- **iter 선정 기준:** CD_occ가 plateau 95%에 도달하는 최소 iter

---

---

# APPENDIX: 전체 선행연구 조사

---

## A. 연구 동향과 방법론 계보

### A.1 장애물 표현 방식의 분류

점군 기반 장애물 표현은 크게 두 갈래로 나뉜다.

**명시적 표현 (Explicit)**
메시(mesh), 점유 격자(occupancy grid), DEM(Digital Elevation Model) 등이 해당한다. 표면 좌표를 직접 저장하므로 시각화에는 우수하나, 내부/외부 판정에 ray-casting 같은 추가 처리가 필요하고 GNC 제약으로 직접 삽입이 불가능하다. 우주 도메인에서는 OSIRIS-REx OLA, Hera PALT 같은 레이저 고도계가 DEM을 생성하지만 지질 해석용이며 실시간 회피 최적화에는 부적합하다.

**음함수 표현 (Implicit)**
$f(\mathbf{x}) = 0$ 의 등위면으로 표면을 정의한다. 부호 하나로 내부/외부를 즉시 판정하고, $\nabla f$ 가 폐쇄형으로 존재하여 GNC 제약에 직접 삽입 가능하다.

---

### A.2 학습 기반 음함수 (Neural Implicit)

**NeRF 계열**
- Pantic et al. (arXiv 2022): NeRF를 장애물 그래디언트 소스로 활용한 반응 계획
- Chen et al. (T-RO 2024, CATNIPS): NeRF → Poisson Point Process 변환, 충돌 확률 정량화

**신경망 SDF 계열**
- Ortiz et al. (RSS 2022, iSDF): MLP가 3D 좌표 → SDF 값을 온라인 학습
- Zhou et al. (T-AES 2025): 비협조 우주 타겟 LiDAR 점군에서 신경망 SDF 학습 + 2층 CBF. **현재 우주 도메인 SOTA**. 단, MLP 파라미터 수만~수십만 개, 폐쇄형 그래디언트 없음.

**공통 한계:** 학습 데이터 의존성, 파라미터 수 수만~수십만, 분석적 폐쇄형 그래디언트 부재.

---

### A.3 해석 기반 음함수 (Analytic Implicit)

#### A.3.1 가우시안 혼합 모델 (GMM)

| 논문 | 저널 | 방법 | 한계 |
|---|---|---|---|
| Chen et al. 2021 | Adv. Space Res. | GMM-APF (회전 타겟) | SOS 불가, 점군 직접 피팅 아님 |
| Cao et al. 2020 | Applied Sciences | GMM-APF + 고정시간 제어 | 동일 |

#### A.3.2 초이차 (Superquadric)

| 논문 | 저널 | 방법 | 한계 |
|---|---|---|---|
| Badawy & McInnes 2008 | JGCD | 초이차 APF + SMC | CAD 파라미터 필요 |
| Park & D'Amico 2024 | SciTech | CNN 2D→3D 초이차 | 점군 직접 피팅 아님 |

#### A.3.3 음함수 다항식 (Implicit Polynomial, IP)

단일 다항식 등위면 $f(\mathbf{x}) = 1$ 로 장애물 표면을 정의한다. 차수 4에서 28개(동차), 차수 6에서 84개(비동차)/28개(동차)의 계수만으로 표현되며 폐쇄형 그래디언트·헤시안이 존재한다.

**피팅 방법론의 두 갈래**
- 대수적 방법(Algebraic): $\sum_i f(\mathbf{x}_i)^2$ 최소화, 선형 시스템, 빠름. 3L 알고리즘이 대표.
- 기하학적 방법(Geometric): 직교 거리 최소화, 비선형 솔버, 정확하지만 느림.

**수치 불안정성과 기존 해결책**
- Ridge Regression: Tasdizen et al. (TIP 2000), 사후 페널티로 조건수 완화
- 제약 IP: Keren & Gotsman (TPAMI 1999), 닫힌 곡면 보장
- 음함수 B-spline (IBS): Rouhani & Sappa (CVPR 2010), 국소 기저로 희소성 확보

세 전략 모두 단항식 기저를 고정한 채 사후에 문제를 해결.

**SOS 확장**
Ahmadi et al. (RSS 2017): SOS 부수준집합으로 두 객체 간 유클리드 거리를 SDP로 계산. 경로계획기 삽입 가능. 단, 정적 메시 가정, 이상치 처리 없음, 우주 도메인 미적용.

**지상 도메인 실시간 확장**
de Sa et al. (ICRA 2024): 이차 IP를 OLS로 피팅 → CBF-QP 직접 투입. 차수 2의 형상 표현 한계, 이상치 처리 없음.

| 논문 | 저널 | 방법 | 주요 기여 |
|---|---|---|---|
| Blane et al. 2000 | TPAMI | 3L (대수적) | IP 피팅 실용화 기준선 |
| Tasdizen et al. 2000 | TIP | Ridge IP | 유클리드 불변 정규화 |
| Keren & Gotsman 1999 | TPAMI | Constrained IP | 위상 보장 |
| Rouhani & Sappa 2010 | CVPR | IBS (B-spline 기저) | 국소 기저, 조건수 개선 |
| Zheng et al. 2010 | TPAMI | Adaptive geometric | 직교 거리 기반 개선 |
| Ahmadi et al. 2017 | RSS | SOS + SDP | GNC 제약 삽입 가능성 |
| de Sa et al. 2024 | ICRA | 이차 IP + CBF-QP | 실시간 지상 적용 |

---

### A.4 우주 도메인 종합

| 논문 | 저널 | 방법 | 목적 |
|---|---|---|---|
| Zhao et al. 2018 | Sensors | ICP + 적응 복셀 | 자세 추적 |
| Renaut et al. 2024 | T-AES | CNN + LiDAR 점군 | 자세 추정 |
| Vela-Rincón et al. 2022 | Acta Astro. | 메시 기반 GNC | GNC 설계 |
| Leomanni et al. 2022 | JGCD | 다면체 LP-MPC | 경로계획 |
| Zhou et al. 2025 | T-AES | 신경망 SDF + CBF | 안전 제어 (SOTA) |

**핵심 공백:** "3D 점군 → 비신경망 해석적 음함수 다항식 피팅 → GNC 제약"의 루프를 비협조 우주 도메인에서 완결한 동료심사 논문은 현재 존재하지 않는다.

---

## B. 본 연구의 방법론 (HERO)

### B.1 Fischer 분해 기반 직교 구면조화 기저

$$f_k(\mathbf{x}) = \sum_{j=0}^{\lfloor k/2 \rfloor} r^{2j} H_{k-2j}(\mathbf{x}), \quad \Delta H_m(\mathbf{x}) = 0$$

- 단위 구 위에서 직교 정규화 → 설계행렬 조건수 구조적 감소
- 사후 Ridge 없이도 피팅 안정성 확보 (단항식 기저 대비 핵심 차이)
- SO(3) 회전 하에서 동일 차수 ℓ 내부에서만 Wigner D-행렬로 변환 → 계수가 시불변

### B.2 Sobolev 가중 Ridge

$$D = \lambda \cdot \text{diag}\left(\frac{[l(l+1)]^s}{[k(k+1)]^s}\right)$$

- l=0 (상수항): 페널티 0 (자동 면제)
- l=2: 페널티 6/42·λ
- l=6: 페널티 λ (최대)
- 물리적 의미: 고차 구면조화(세밀한 형상)를 더 강하게 억제 → 부분 관측 시 과적합 방지

### B.3 RANSAC 2단 스코어링

```
Stage 1 (Raw):
  idx ← randsample(1:N, k·nT, ω 기반)
  β_tmp ← regressionFourthOrder(coord[idx], Φ, λ_Sobolev)
  S_j = mean(exp(-r²/2σ²))

Stage 2 (Local refinement) — S_j ≥ S_best × τ 일 때:
  [δ, β_loc] ← DisplacementLocal(inlier 점, Φ, ∇Φ)
  S_jLoc = mean(exp(-r_loc²/2σ²))
  ω ← (S_sum·ω + S_jLoc·ω_jLoc) / (S_sum + S_jLoc)
```

### B.4 기여 요약 (C1~C3)

| 기여 | 내용 | 선행 연구 대비 |
|---|---|---|
| C1 | 짝수차 동차다항식 → 비협조 우주 객체 최초 적용 | 기존 우주 연구: 자세 추정, GNC 연결 없음 |
| C2 | Fischer 기저 → κ 구조적 감소 + Sobolev Ridge 정당화 | 기존 IP: 단항식 기저 + 사후 Ridge |
| C3 | RANSAC + 적응 재가중치 → IP 분야 최초 명시적 이상치 거부 | 기존: 깨끗한 메시 가정 |
| C4 | KNN 그래프 + 스펙트럴 클러스터링 전처리 → 수동 레이블 없이 다물체 구조 분리 후 HERO 적용 | 기존: 단일 객체 가정 또는 수동 분리 |

**C4 포지셔닝:** 클러스터링 자체(기존 기법)가 기여가 아니라, **HERO 파이프라인이 전처리와 결합하여 복잡 다물체 우주 구조물로 자연스럽게 확장됨**이 기여. 졸논 구현을 이식하여 적용.

---

## C. 선행연구 상세 조사 (점군 피팅·표현 방식별)

### C.1 음함수 다항식 (Implicit Polynomial) 계열

#### Blane et al. 2000 — 3L 알고리즘

| 항목 | 내용 |
|---|---|
| 저자 | M. M. Blane, T. Lei, H. Çivi, D. B. Cooper |
| 제목 | The 3L Algorithm for Fitting Implicit Polynomial Curves and Surfaces to Data |
| 게재 | IEEE TPAMI, vol. 22, no. 3, pp. 298–313 |
| DOI | 10.1109/34.841760 |

3D 데이터: Princeton Shape Benchmark 3D 메시 샘플링 점군. 경계면 내·외부 두 오프셋 수준집합으로 IP 피팅을 선형 시스템으로 정식화. IP 피팅 실용화 기준선.

#### Tasdizen et al. 2000 — Ridge IP

| 항목 | 내용 |
|---|---|
| 저자 | T. Tasdizen, J.-P. Tarel, D. B. Cooper |
| 제목 | Improving the Stability of Algebraic Curves for Applications |
| 게재 | IEEE TIP, vol. 9, no. 3, pp. 405–416 |
| DOI | 10.1109/83.826778 |

유클리드 불변 3D Ridge 행렬 명시적 유도. 사후 페널티로 조건수 완화. 단순 단위 행렬 Ridge 대비 좌표계 방향 무관한 일관된 정규화.

#### Keren & Gotsman 1999 — 제약 IP

| 항목 | 내용 |
|---|---|
| 저자 | D. Keren, C. Gotsman |
| 제목 | Fitting Curves and Surfaces With Constrained Implicit Polynomials |
| 게재 | IEEE TPAMI, vol. 21, no. 1, pp. 31–41 |
| DOI | 10.1109/34.745731 |

영집합이 star-shaped 또는 유계가 되도록 위상학적 제약 부과. 닫힌 곡면 보장. 위상 불일치(phantom sheets) 방지.

#### Rouhani & Sappa 2010 — IBS

| 항목 | 내용 |
|---|---|
| 저자 | M. Rouhani, A. D. Sappa |
| 제목 | Implicit B-Spline Fitting Using the 3L Algorithm |
| 게재 | IEEE CVPR 2010 |

전역 단항식 기저 → 격자 기반 B-spline 국소 기저로 교체. 설계행렬 구조적 희소화, 조건수 개선. tension 정규화.

#### Zheng et al. 2010 — Adaptive Geometric

| 항목 | 내용 |
|---|---|
| 저자 | J. Zheng, J. Takamatsu, K. Ikeuchi |
| 제목 | An Adaptive and Stable Method for Fitting Implicit Polynomial Curves and Surfaces |
| 게재 | IEEE TPAMI, vol. 32, no. 3, pp. 561–568 |
| DOI | 10.1109/TPAMI.2009.189 |

직교 거리 최소화 기반 기하적 방법. Levenberg–Marquardt 수렴 안정성 개선. 계산 비용 큼.

#### Ahmadi et al. 2017 — SOS + SDP ★ 본 연구 직접 선조

| 항목 | 내용 |
|---|---|
| 저자 | A. A. Ahmadi, G. Hall, A. Makadia, V. Sindhwani |
| 제목 | Geometry of 3D Environments and Sum of Squares Polynomials |
| 게재 | RSS XIII |
| DOI | 10.15607/RSS.2017.XIII.071 |

SOS 부수준집합으로 두 객체 간 유클리드 거리를 SDP로 계산. 경로계획기 삽입 가능한 폐쇄형 거리 oracle. 한계: 정적 메시 가정, 이상치 없음, 우주 미적용.

#### de Sa et al. 2024 — 이차 IP + CBF-QP ★ 본 연구 직접 비교군

| 항목 | 내용 |
|---|---|
| 저자 | M. de Sa, P. Kotaru, K. Sreenath |
| 제목 | Point Cloud-Based Control Barrier Function Regression for Safe and Efficient Vision-Based Control |
| 게재 | IEEE ICRA 2024 |
| DOI | 10.1109/ICRA57147.2024.10610647 |

매 주기 이차 IP OLS 피팅 → CBF-QP. SLAM 없이 실시간 안전 보장. 한계: 차수 2, 이상치 처리 없음, 우주 미적용.

---

### C.2 신경망 SDF / NeRF 계열

#### Ortiz et al. 2022 — iSDF

| 항목 | 내용 |
|---|---|
| 저자 | J. Ortiz et al. |
| 제목 | iSDF: Real-Time Neural Signed Distance Fields for Robot Perception |
| 게재 | RSS 2022 |
| arXiv | 2204.02296 |

MLP 온라인 학습: 3D 좌표 → SDF. Voxel Grid 대비 적응적 해상도, 컴팩트.

#### Chen et al. 2024 — CATNIPS

| 항목 | 내용 |
|---|---|
| 저자 | T. Chen et al. |
| 제목 | CATNIPS: Collision Avoidance Through Neural Implicit Probabilistic Scenes |
| 게재 | IEEE T-RO, vol. 40, pp. 2712–2733 |
| DOI | 10.1109/TRO.2024.3387428 |

NeRF → Poisson Point Process 변환, 충돌 확률 정량화.

#### Zhou et al. 2025 — 신경망 SDF + CBF (우주 도메인 SOTA) 🛰

| 항목 | 내용 |
|---|---|
| 저자 | Y. Zhou, Y. Shi, H. Mao et al. |
| 제목 | Spacecraft Safe Robust Control Using Implicit Neural Representation for Geometrically Complex Targets in Proximity Operations |
| 게재 | IEEE T-AES 2025 |
| arXiv | 2507.13672 |

비협조 우주 타겟 LiDAR 점군에서 직접 SDF 학습. 2층 계층적 CBF. 한계: MLP 파라미터 수만~수십만 개, 폐쇄형 그래디언트 없음.

---

### C.3 기하 프리미티브 계열 (초이차·GMM)

#### Badawy & McInnes 2008 🛰

| 항목 | 내용 |
|---|---|
| 저자 | A. Badawy, C. R. McInnes |
| 제목 | On-Orbit Assembly Using Superquadric Potential Fields |
| 게재 | JGCD, vol. 31, no. 1, pp. 30–43 |
| DOI | 10.2514/1.28865 |

초이차 포텐셜 함수 기반 우주 근접 운용 APF 최초 제안. CAD 파라미터 필요.

#### Chen et al. 2021 🛰

| 항목 | 내용 |
|---|---|
| 저자 | X. Chen et al. |
| 제목 | Obstacle Avoidance for Non-Cooperative Target Spacecraft with Gaussian Mixture Model |
| 게재 | Adv. Space Res., vol. 68, no. 10, pp. 4217–4233 |
| DOI | 10.1016/j.asr.2021.08.009 |

GMM 형태 척력 포텐셜. 회전 타겟 대응. SOS 불가.

#### Cao et al. 2020 🛰

| 항목 | 내용 |
|---|---|
| 저자 | C. Cao et al. |
| 제목 | Obstacle Avoidance for Spacecraft in Close Proximity Operations Using GMM and Fixed-Time Control |
| 게재 | Applied Sciences, vol. 10, no. 17, p. 5986 |
| DOI | 10.3390/app10175986 |

GMM-APF + 고정시간 수렴 제어. 임무 시간 제약 만족.

#### Park & D'Amico 2024 🛰

| 항목 | 내용 |
|---|---|
| 저자 | T. H. Park, S. D'Amico |
| 제목 | Rapid Abstraction of Spacecraft 3D Structure from Single 2D Image |
| 게재 | AIAA SciTech 2024 |
| DOI | 10.2514/6.2024-0963 |

CNN 단일 2D 이미지 → 3D 초이차 어셈블리 파라미터 예측. 점군 직접 피팅 아님.

#### Liu et al. 2022 — 확률적 초이차

| 항목 | 내용 |
|---|---|
| 저자 | L. Liu et al. |
| 제목 | Robust and Accurate Superquadric Recovery: A Probabilistic Approach |
| 게재 | CVPR 2022 |

GMM + EM 프레임워크로 초이차 피팅을 MLE 문제로 정식화. 이상치 오염 60%까지 강건.

---

### C.4 우주 도메인 — 자세 추정 (장애물 표현 목적 아님)

| 논문 | 저널 | 방법 | 용도 |
|---|---|---|---|
| Zhao et al. 2018 | Sensors | ICP + 적응 복셀 | 자세 추적 |
| Renaut et al. 2024 | T-AES | CNN + LiDAR | 자세 추정 (대칭 처리) |
| Vela-Rincón et al. 2022 | Acta Astro. | 메시 기반 GNC | GNC 설계 |
| Leomanni et al. 2022 | JGCD | 다면체 KOZ + LP-MPC | 경로계획 |

---

### C.5 지상/드론 — 점군 기반 장애물 처리

#### Dai et al. 2024 — Sailing CBF

| 항목 | 내용 |
|---|---|
| 저자 | B. Dai et al. |
| 제목 | Sailing Through Point Clouds: Safe Navigation Using Point Cloud-Based Control Barrier Functions |
| 게재 | IEEE RA-L, vol. 9, no. 9, pp. 7731–7738 |
| DOI | 10.1109/LRA.2024.3431870 |

점군에서 직접 스케일링 인수 추출 → CBF-QP. 피팅 단계 없이 직접 활용.

#### Wen et al. 2024 — Implicit Swept Volume SDF

| 항목 | 내용 |
|---|---|
| 저자 | W. Wen et al. |
| 제목 | Implicit Swept Volume SDF |
| 게재 | ACM TOG (SIGGRAPH), vol. 43, no. 4 |
| DOI | 10.1145/3658181 |

임의 형상 로봇의 swept volume을 신경망 SDF로 연속 표현.

---

## D. 비교 종합

### D.1 표현 방식별 핵심 트레이드오프

| 표현 방식 | 대표 논문 | 파라미터 수 | 해석적 그래디언트 | SDP 호환 | 이상치 강건 | 우주 도메인 |
|---|---|---|---|---|---|---|
| Monomial IP (차수 6, 동차) | Blane 2000 | 28 | ✓ | ✓ | △ | ✗ |
| 이차 IP (CBF 회귀) | de Sa 2024 | 10 | ✓ | ✓ | ✗ | ✗ |
| 초이차 | Park 2024 | 5–15 | ✓ | △ | △ | 🛰 일부 |
| GMM | Chen 2021 | ~90 | ✓ | ✗ | △ | 🛰 일부 |
| 신경망 SDF | Zhou 2025 | 수만~수십만 | △ (자동미분) | ✗ | △ | 🛰 Zhou 2025 |
| **HERO (본 연구)** | **정우진 2026** | **28 (동차 6차)** | **✓ (폐쇄형)** | **✓** | **✓ (RANSAC)** | **🛰 최초** |

### D.2 핵심 리서치 갭

우주 도메인 논문(Zhao 2018, Renaut 2024, Vela-Rincón 2022)은 모두 3D LiDAR 점군을 사용하나 목적이 자세 추정이며, 해석적 장애물 표면 모델을 피팅하여 GNC 제약으로 연결하지 않는다. Zhou 2025가 유일하게 이 루프를 구성하나 신경망 SDF 기반이다.

**"3D 점군 → 비신경망 음함수 다항식 피팅 → GNC 제약"의 루프를 비협조 우주 도메인에서 완결한 동료심사 논문은 현재까지 존재하지 않는다.**

---

## E. 참고문헌

### IP 피팅 계보

[1] M. M. Blane, T. Lei, H. Çivi, D. B. Cooper, "The 3L Algorithm for Fitting Implicit Polynomial Curves and Surfaces to Data," *IEEE TPAMI*, vol. 22, no. 3, pp. 298–313, 2000.

[2] T. Tasdizen, J.-P. Tarel, D. B. Cooper, "Improving the Stability of Algebraic Curves for Applications," *IEEE TIP*, vol. 9, no. 3, pp. 405–416, 2000.

[3] D. Keren, C. Gotsman, "Fitting Curves and Surfaces With Constrained Implicit Polynomials," *IEEE TPAMI*, vol. 21, no. 1, pp. 31–41, 1999.

[4] M. Rouhani, A. D. Sappa, "Implicit B-Spline Fitting Using the 3L Algorithm," *CVPR*, 2010.

[5] J. Zheng, J. Takamatsu, K. Ikeuchi, "An Adaptive and Stable Method for Fitting Implicit Polynomial Curves and Surfaces," *IEEE TPAMI*, vol. 32, no. 3, pp. 561–568, 2010.

[6] A. A. Ahmadi, G. Hall, A. Makadia, V. Sindhwani, "Geometry of 3D Environments and Sum of Squares Polynomials," *RSS XIII*, 2017.

[7] M. de Sa, P. Kotaru, K. Sreenath, "Point Cloud-Based CBF Regression for Safe and Efficient Vision-Based Control," *ICRA*, 2024.

### 신경망 음함수

[8] J. Ortiz et al., "iSDF: Real-Time Neural Signed Distance Fields for Robot Perception," *RSS*, 2022.

[9] T. Chen et al., "CATNIPS: Collision Avoidance Through Neural Implicit Probabilistic Scenes," *IEEE T-RO*, vol. 40, 2024.

[10] Y. Zhou et al., "Spacecraft Safe Robust Control Using Implicit Neural Representation for Geometrically Complex Targets," *IEEE T-AES*, 2025.

### 우주 도메인

[11] A. Badawy, C. R. McInnes, "On-Orbit Assembly Using Superquadric Potential Fields," *JGCD*, vol. 31, no. 1, 2008.

[12] X. Chen et al., "Obstacle Avoidance for Non-Cooperative Target Spacecraft with GMM," *Adv. Space Res.*, vol. 68, no. 10, 2021.

[13] C. Cao et al., "Obstacle Avoidance Using GMM and Fixed-Time Control," *Applied Sciences*, vol. 10, no. 17, 2020.

[14] M. Leomanni et al., "Explicit MPC for Low-Thrust Spacecraft Proximity Operations," *JGCD*, vol. 45, no. 5, 2022.

[15] T. H. Park, S. D'Amico, "Rapid Abstraction of Spacecraft 3D Structure from Single 2D Image," *AIAA SciTech*, 2024.

[16] A. Vela-Rincón et al., "GNC Design for Proximity Operations with Tumbling Target," *Acta Astro.*, vol. 196, 2022.

[17] G. Zhao et al., "LiDAR-Based Non-Cooperative Tumbling Spacecraft Pose Tracking," *Sensors*, vol. 18, no. 10, 2018.

[18] L. Renaut et al., "CNN-based Pose Estimation of Non-Cooperative Spacecraft from LiDAR," *IEEE T-AES*, 2024.

### 3D 점군 피팅 및 표면 재구성

[19] B. Dai et al., "Sailing Through Point Clouds: Safe Navigation Using PCB-CBF," *IEEE RA-L*, vol. 9, no. 9, 2024.

[20] W. Wen et al., "Implicit Swept Volume SDF," *ACM TOG (SIGGRAPH)*, vol. 43, no. 4, 2024.

[21] L. Liu et al., "Robust and Accurate Superquadric Recovery: A Probabilistic Approach," *CVPR*, 2022.

---

*최종 업데이트: 2026년 5월*
