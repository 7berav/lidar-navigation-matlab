# 논문 확장 전략 종합 요약

**대상:** Even-Degree Homogeneous Polynomial Modeling for 3D Obstacle Representation (정우진)
**목적:** T-AES Navigation / ASR 투고를 위한 연구 범위·방법론·지표 정리
**작성일:** 2026년 5월

---

## 0. 핵심 결론 (한 줄)

"3D LiDAR 점군 → Fischer 분해 기반 짝수차 동차다항식 음함수 피팅 → GNC 제약 준비"의 파이프라인을 비협조 우주 도메인에서 최초로 완결한 논문으로, **피팅 방법론 단독으로 1편, 제어 통합으로 2편**의 2단 출판 전략을 권고함.

---

## 1. 연구 동향과 방법론 계보

### 1.1 장애물 표현 방식의 분류

점군 기반 장애물 표현은 크게 두 갈래로 나뉜다.

**명시적 표현 (Explicit)**
메시(mesh), 점유 격자(occupancy grid), DEM(Digital Elevation Model) 등이 해당한다. 표면 좌표를 직접 저장하므로 시각화에는 우수하나, 내부/외부 판정에 ray-casting 같은 추가 처리가 필요하고 GNC 제약으로 직접 삽입이 불가능하다. 우주 도메인에서는 OSIRIS-REx OLA, Hera PALT 같은 레이저 고도계가 DEM을 생성하지만 지질 해석용이며 실시간 회피 최적화에는 부적합하다.

**음함수 표현 (Implicit)**
$f(\mathbf{x}) = 0$ 의 등위면으로 표면을 정의한다. 부호 하나로 내부/외부를 즉시 판정하고, $\nabla f$ 가 폐쇄형으로 존재하여 GNC 제약에 직접 삽입 가능하다. 음함수 방법론은 다시 두 가지로 나뉜다.

---

### 1.2 학습 기반 음함수 (Neural Implicit)

**NeRF 계열**
Pantic et al. (arXiv 2022)은 NeRF를 장애물 그래디언트 소스로 활용한 반응 계획을 제안하였다. Chen et al. (T-RO 2024, CATNIPS)은 NeRF를 Poisson Point Process로 변환하여 충돌 확률을 정량화하고 드론 궤적 계획에 적용하였다.

**신경망 SDF 계열**
Ortiz et al. (RSS 2022, iSDF)은 깊이 카메라 점군 스트림으로부터 MLP가 3D 좌표 → SDF 값을 온라인 학습하는 구조를 제안하였다. Zhou et al. (T-AES 2025)은 이를 우주 도메인으로 확장하여 비협조 타겟 LiDAR 점군에서 신경망 SDF를 학습하고 2층 CBF 안전 필터를 구성하였다. 현재 비협조 우주 근접 운용 분야의 SOTA이나, MLP 파라미터가 수만 개에 달하여 방사선 강화(radiation-hardened) 비행 컴퓨터에 부담이 크고 SDP/QP에 직접 삽입 가능한 폐쇄형 그래디언트가 없다.

**공통 한계:** 학습 데이터 의존성, 파라미터 수 수만~수십만, 분석적 폐쇄형 그래디언트 부재.

---

### 1.3 해석 기반 음함수 (Analytic Implicit)

#### 1.3.1 가우시안 혼합 모델 (GMM)

Chen et al. (ASR 2021)은 비협조 타겟의 3D 형상에 GMM을 적합하고 각 가우시안 커널로 척력 포텐셜을 구성하였다. Cao et al. (Applied Sciences 2020)은 GMM-APF에 고정시간(fixed-time) 수렴 제어를 결합하여 임무 시간 제약을 만족하는 근접 운용을 시연하였다. 회전 타겟에 부분적으로 대응할 수 있으나, 다항식 형태가 아니므로 SOS 거리 인증이 불가능하다.

| 논문 | 저널 | 방법 | 한계 |
|---|---|---|---|
| Chen et al. 2021 | Adv. Space Res. | GMM-APF (회전 타겟) | SOS 불가, 점군 직접 피팅 아님 |
| Cao et al. 2020 | Applied Sciences | GMM-APF + 고정시간 제어 | 동일 |

#### 1.3.2 초이차 (Superquadric)

Badawy & McInnes (JGCD 2008)는 초이차 포텐셜 함수를 이용한 우주 근접 운용 APF를 최초로 제안하였다. Park & D'Amico (SciTech 2024)는 단일 2D 이미지에서 CNN으로 3D 초이차 어셈블리 파라미터를 예측하였다. 해석적이고 컴팩트하나, 형상 파라미터를 CAD 모델 또는 학습 데이터로부터 설정해야 하므로 사전 정보 없이 LiDAR 점군에서 직접 피팅하는 것이 어렵다.

| 논문 | 저널 | 방법 | 한계 |
|---|---|---|---|
| Badawy & McInnes 2008 | JGCD | 초이차 APF + SMC | CAD 파라미터 필요 |
| Park & D'Amico 2024 | SciTech | CNN 2D→3D 초이차 | 점군 직접 피팅 아님 |

#### 1.3.3 음함수 다항식 (Implicit Polynomial, IP)

단일 다항식 등위면 $f(\mathbf{x}) = 1$ 로 장애물 표면을 정의한다. 차수 4에서 35개, 차수 6에서 84개의 계수만으로 표현되며 폐쇄형 그래디언트·헤시안이 존재한다. Sum-of-Squares(SOS) 프레임워크와 직접 결합 가능하여 거리 인증이 가능하다.

**피팅 방법론의 두 갈래**

대수적 방법(Algebraic)은 $\sum_i f(\mathbf{x}_i)^2$ 를 최소화하는 선형 시스템을 풀며, 빠르고 구현이 단순하지만 실제 유클리드 거리의 의미가 없어 정확도가 한정된다. 3L 알고리즘(Blane et al. TPAMI 2000)이 대표적이며, 경계면 내·외부에 두 오프셋 수준집합을 구성하여 부호 전환 제약을 선형 시스템으로 정식화한다.

기하학적 방법(Geometric)은 $\sum_i d(\mathbf{x}_i, \mathcal{S})^2$ 즉 직교 거리를 최소화하며 Levenberg–Marquardt 등 반복 비선형 솔버를 사용한다. 정확도가 높으나 계산 비용이 크다.

**수치 불안정성과 기존 해결책**

고차 단항식 기저는 설계 행렬 열 간 강한 상관으로 조건수가 극단적으로 커진다. 기존 해결 전략은 세 가지다.

- Ridge Regression: Tasdizen et al. (TIP 2000)이 유클리드 불변 3D Ridge 행렬을 유도하여 사후(post-hoc) 페널티 항으로 조건수를 완화하였다.
- 제약 IP: Keren & Gotsman (TPAMI 1999)이 영집합의 위상 형태를 제약하여 닫힌 곡면을 보장하였다.
- 음함수 B-spline(IBS): Rouhani & Sappa (CVPR 2010)이 전역 단항식 기저를 격자 기반 B-spline 국소 기저로 교체하여 구조적 희소성과 조건수를 개선하였다.

세 전략 모두 단항식 기저를 고정한 채 사후에 문제를 해결한다는 공통점이 있다.

**SOS 확장**

Ahmadi et al. (RSS 2017)은 SOS 부수준집합을 이용해 두 객체 간 유클리드 거리를 SDP로 계산하는 거리 오라클을 제안하였다. 경로계획기 삽입 가능성을 시연하였으나, 정적 메시 가정으로 텀블링 타겟에 미적용이고 이상치 처리가 없다.

**지상 도메인 실시간 확장**

de Sa et al. (ICRA 2024)은 매 제어 주기마다 RealSense 깊이 카메라 점군에 이차 IP를 OLS로 피팅하여 CBF-QP에 직접 투입하는 방식을 제안하였다. SLAM 없이 실시간 안전 보장을 달성하였으나, 차수 2 이차식의 형상 표현 한계와 이상치 명시 처리 부재가 한계다.

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

### 1.4 우주 도메인 종합

우주 근접 운용 분야에서 3D 점군을 직접 사용하는 연구들은 대부분 **자세 추정(pose estimation)** 목적이다.

| 논문 | 저널 | 방법 | 목적 |
|---|---|---|---|
| Zhao et al. 2018 | Sensors | ICP + 적응 복셀 | 자세 추적 |
| Renaut et al. 2024 | T-AES | CNN + LiDAR 점군 | 자세 추정 |
| Vela-Rincón et al. 2022 | Acta Astro. | 메시 기반 GNC | GNC 설계 |
| Leomanni et al. 2022 | JGCD | 다면체 LP-MPC | 경로계획 |
| Zhou et al. 2025 | T-AES | 신경망 SDF + CBF | 안전 제어 (SOTA) |

**핵심 공백:** "3D 점군 → 비신경망 해석적 음함수 다항식 피팅 → GNC 제약"의 루프를 비협조 우주 도메인에서 완결한 동료심사 논문은 현재 존재하지 않는다.

---

## 2. 본 연구의 방법론과 기여

### 2.1 문제 설정

비협조 우주 타겟의 3D LiDAR 점군을 입력으로 받아, 노이즈와 이상치에 강건하게 짝수차 동차다항식 음함수 표면을 추정하는 파이프라인을 제안한다. 대상은 중심대칭(centrally symmetric) 물체로 한정한다. 짝수차 동차다항식은 $f(-\mathbf{x}) = f(\mathbf{x})$ 를 수학적으로 만족하며, 이는 구조·질량 균형 요건으로 중심대칭 설계가 보편적인 우주선에 자연스럽게 부합한다.

### 2.2 핵심 방법론

**Fischer 분해 기반 직교 구면조화 기저**

기존 단항식 $x^a y^b z^c$ 기저 대신, Fischer 분해를 통해 동차다항식을 방사 인수(radial factor)와 구면조화(solid spherical harmonics)의 합으로 분해한 직교 기저를 사용한다.

$$f_k(\mathbf{x}) = \sum_{j=0}^{\lfloor k/2 \rfloor} r^{2j} H_{k-2j}(\mathbf{x}), \quad \Delta H_m(\mathbf{x}) = 0$$

이 기저는 단위 구 위에서 직교 정규화되므로 설계 행렬의 조건수가 구조적으로 낮아진다. 사후 Ridge 페널티 없이 피팅 안정성이 확보되는 것이 단항식 기저 대비 핵심 차이다. 추가로, 구면조화는 SO(3) 회전 하에서 동일 차수 $\ell$ 내부에서만 Wigner D-행렬로 섞인다. 이로 인해 타겟 기체 프레임의 다항식 계수는 타겟 텀블링과 무관하게 시불변이며, 관성 프레임으로의 변환이 해석적으로 전파된다.

**RANSAC 기반 강건 추정 (페널티 반영 최적화)**

우주 LiDAR 점군은 태양광 정반사, 표면 반사율 변화, 부분 가림으로 인한 이상치를 포함한다. 매 반복에서 무작위 서브샘플로 피팅 가설을 생성하고, 잔차에 기반한 적응 재가중치(adaptive residual reweighting)로 이상치를 명시적으로 거부한다. 이는 기존 IP 피팅 연구에서 사실상 최초로 RANSAC을 결합한 사례다.

목적함수는 가중 잔차 최소화 형태이다.

$$\min_{\boldsymbol{\beta}} \sum_{i} w_i \left( f(\mathbf{x}_i; \boldsymbol{\beta}) - 1 \right)^2, \quad w_i = \text{score}(\text{residual}_i)$$

가중치 $w_i$ 는 잔차 크기에 따라 적응적으로 감소하므로, 이상치의 영향이 반복 과정에서 자동으로 줄어든다.

**다물체 스펙트럴 클러스터링**

ISS 같은 복잡 타겟은 모듈, 솔라패널, 트러스가 위상학적으로 분리된 구조를 가지므로 단일 다항식으로 표현이 불가능하다. 점간 거리와 표면 법선 유사도로 근접 그래프를 구성하고 정규화 라플라시안의 스펙트럴 클러스터링으로 구성 요소를 분리한 후, 각 클러스터에 독립적인 HERO 모델을 피팅한다.

### 2.3 기여 요약

| 기여 | 내용 | 선행 연구 대비 차이 |
|---|---|---|
| C1: 표현 | 짝수차 동차다항식 음함수 등위면 → 비협조 우주 객체 최초 적용 | 기존 우주 연구: 자세 추정 목적, GNC 연결 없음 |
| C2: 기저 | Fischer 분해 직교 구면조화 → 구조적 조건수 개선 + R(t) 호환 | 기존 IP 연구: 단항식 기저 + 사후 Ridge |
| C3: 강건성 | RANSAC + 적응 재가중치 → IP 분야 최초 명시적 이상치 거부 | 기존: 깨끗한 메시 가정 또는 시간 평균에 의존 |
| C4: 다물체 | 스펙트럴 클러스터링 + 컴포넌트별 피팅 → 수동 레이블 없이 | 기존: 단일 객체 가정 또는 수동 분리 |

---

## 3. 시나리오 및 검증 계획

### 3.1 투고 전략 (2단 로켓)

**1편 — 피팅 방법론 단독 (현재 투고 가능)**
- 타겟 저널: T-AES Navigation 영역 또는 Advances in Space Research
- 내용: C1~C4 기여 + 피팅 성능 검증
- 제어 실험 불필요. 마지막 절에 "폐루프 GNC 통합은 향후 연구"로 명시

**2편 — 제어 통합 (1편 게재 후 6~12개월)**
- 타겟 저널: T-AES Autonomous Systems 또는 JGCD
- 내용: 1편의 IP 모델을 CBF-QP 또는 MPC 제약으로 삽입한 폐루프 시뮬레이션
- Zhou 2025와 동일 시나리오에서 정면 비교

### 3.2 검증 시나리오 후보

#### 시나리오 A — 피팅 성능 단독 검증 (1편 핵심)

| 항목 | 내용 |
|---|---|
| 타겟 | CubeSat 형상, ISS 형상 (단일 + 다물체) |
| 입력 | 시뮬레이션 LiDAR 점군 (SISPO 또는 Gazebo 기반) |
| 노이즈 | 이상치 오염 10 / 20 / 30 / 40% 몬테카를로 |
| 비교군 | 신경망 SDF (Zhou 2025 재구현), 이차 IP OLS (de Sa 2024), GMM |

#### 시나리오 B — 단순 GNC 통합 검증 (1편 선택 또는 2편)

비복잡한 기존 방법론을 가져다 쓰는 것이 핵심이다. 제어 기법 자체가 기여가 아니라 "IP 모델이 기존 제어 프레임워크에 플러그인 가능함을 시연"하는 것으로 포지셔닝한다.

**B-1. APF (Artificial Potential Field)**
IP 그래디언트 $\nabla f$ 를 척력 방향으로 직접 사용하는 가장 단순한 방법이다. CW 동역학 위에서 접근 궤적을 생성하며, 구현이 최소화되어 피팅 방법론에 집중할 수 있다.
- 참고: Badawy & McInnes (JGCD 2008), Chen et al. (ASR 2021), Cao et al. (Applied Sciences 2020)

**B-2. CBF-QP**
IP 등위면 $f(\mathbf{x}) \leq 1$ 을 CBF 제약으로 삽입하고 QP를 푸는 방식이다. de Sa et al. (ICRA 2024)의 이차 IP 버전과 직접 비교가 가능하며, 동일 구조에서 차수만 높인 효과를 보일 수 있다. CW 동역학 기반 선형 제어기에 안전 필터로 추가하는 형태로 최소 구현 가능하다.
- 참고: de Sa et al. (ICRA 2024), Zhou et al. (T-AES 2025)

**B-3. 단순 MPC (CW 기반 선형 MPC)**
CW 방정식을 예측 모델로 하고 IP 기반 볼록 근사 제약을 추가하는 선형 MPC이다. SOS 거리를 제약으로 삽입하려면 Ahmadi et al. (RSS 2017)의 SDP 거리 오라클을 그대로 활용할 수 있다.
- 참고: Ahmadi et al. (RSS 2017), Leomanni et al. (JGCD 2022)

### 3.3 유사 연구 — GNC 통합 선례

| 논문 | 저널 | 장애물 표현 | 제어 방법 | 시나리오 |
|---|---|---|---|---|
| Zhou et al. 2025 | T-AES | 신경망 SDF | 2층 CBF | 비협조 위성 근접 |
| Chen et al. 2021 | ASR | GMM | APF | 텀블링 타겟 회피 |
| Cao et al. 2020 | Applied Sci. | GMM | APF + 고정시간 | 텀블링 근접 운용 |
| Leomanni et al. 2022 | JGCD | 다면체 KOZ | LP-MPC | 저추력 근접 |
| Badawy & McInnes 2008 | JGCD | 초이차 | APF + SMC | 궤도상 조립 |
| de Sa et al. 2024 | ICRA | 이차 IP | CBF-QP | 지상/드론 |
| Ahmadi et al. 2017 | RSS | SOS IP | SDP 거리 | 매니퓰레이터/일반 |

---

## 4. 평가 지표 후보

### 4.1 피팅 성능 지표 (1편 핵심)

| 지표 | 설명 | 용도 |
|---|---|---|
| **MSE (잔차 기반)** | $\frac{1}{N}\sum_i (f(\mathbf{x}_i) - 1)^2$ | 피팅 품질 자체 측정. 대수적 잔차이므로 실제 거리 의미 없음. 빠른 계산. |
| **Chamfer Distance** | 두 점군 간 상호 최근접 거리 평균 | 표면 복원 품질 측정. 실제 기하학적 의미 있음. 재구성 평가에 적합. |
| **Hausdorff Distance** | 두 점군 간 최대 최근접 거리 | 최악 오차 측정. Chamfer의 보완 지표. |
| **SDF 오차 (%)** | IP 등위면과 실제 표면 간 거리 / 특성 크기 | 제어 제약으로 쓸 때의 실용적 정확도. Zhou 2025와 직접 비교 가능. |
| **조건수 $\kappa$** | 설계 행렬의 조건수 | Fischer 기저 vs 단항식 기저 수치 안정성 비교. |
| **계산 시간 (ms)** | 피팅 파이프라인 총 소요 시간 | 실시간 온보드 적용 가능성 검증. |
| **메모리 풋프린트 (KB)** | 계수 저장 크기 | 신경망 SDF 대비 압축성 비교. |

**추천 조합:** Chamfer Distance (복원 품질) + SDF 오차 % (GNC 관련성) + 조건수 (Fischer 기저 정당화) + 계산 시간 (실용성)

MSE와 Chamfer를 둘 다 보고할 경우: MSE는 최적화 목적함수와 직접 연결되므로 학습 성능 관점, Chamfer는 실제 기하 복원 품질 관점으로 상보적이다.

### 4.2 강건성 지표 (몬테카를로)

| 지표 | 설명 |
|---|---|
| **이상치 오염률별 SDF 오차 변화** | 오염 10/20/30/40%에서 오차가 얼마나 증가하는가 |
| **수렴 성공률 (%)** | N회 몬테카를로 중 허용 오차 내 수렴 비율 |
| **RANSAC 반복 횟수** | 수렴까지 필요한 평균 반복 수 |

### 4.3 다물체 분리 지표

| 지표 | 설명 |
|---|---|
| **Precision / Recall** | 클러스터 레이블 정확도 |
| **Adjusted Rand Index (ARI)** | 클러스터링 품질 종합 지표 |
| **컴포넌트별 SDF 오차** | 분리 후 각 컴포넌트 피팅 품질 |

### 4.4 GNC 통합 지표 (2편 또는 1편 선택 추가)

| 지표 | 설명 |
|---|---|
| **최소 안전 거리 (m)** | 접근 궤적에서 표면까지 최소 거리 |
| **충돌 여부 / 충돌률 (%)** | N회 몬테카를로 중 충돌 발생 비율 |
| **CBF 제약 위반 횟수** | 폐루프 시뮬레이션 중 CBF 조건 위반 횟수 |
| **Δv 소비량 (m/s)** | 회피 기동 연료 비용 |
| **목적함수 수렴값** | MPC 비용 함수 최종값 비교 |
| **QP/SDP 계산 시간 (ms)** | 실시간 제어 주기 내 해 가능 여부 |

---

## 5. 참고문헌

### IP 피팅 계보

[1] M. M. Blane, T. Lei, H. Çivi, D. B. Cooper, "The 3L Algorithm for Fitting Implicit Polynomial Curves and Surfaces to Data," *IEEE Trans. Pattern Analysis and Machine Intelligence*, vol. 22, no. 3, pp. 298–313, 2000. DOI: 10.1109/34.841760.

[2] T. Tasdizen, J.-P. Tarel, D. B. Cooper, "Improving the Stability of Algebraic Curves for Applications," *IEEE Trans. Image Processing*, vol. 9, no. 3, pp. 405–416, 2000. DOI: 10.1109/83.826778.

[3] D. Keren, C. Gotsman, "Fitting Curves and Surfaces With Constrained Implicit Polynomials," *IEEE Trans. Pattern Analysis and Machine Intelligence*, vol. 21, no. 1, pp. 31–41, 1999. DOI: 10.1109/34.745731.

[4] M. Rouhani, A. D. Sappa, "Implicit B-Spline Fitting Using the 3L Algorithm," in *Proc. IEEE CVPR*, 2010.

[5] J. Zheng, J. Takamatsu, K. Ikeuchi, "An Adaptive and Stable Method for Fitting Implicit Polynomial Curves and Surfaces," *IEEE Trans. Pattern Analysis and Machine Intelligence*, vol. 32, no. 3, pp. 561–568, 2010. DOI: 10.1109/TPAMI.2009.189.

[6] A. A. Ahmadi, G. Hall, A. Makadia, V. Sindhwani, "Geometry of 3D Environments and Sum of Squares Polynomials," in *Proc. RSS XIII*, 2017. DOI: 10.15607/RSS.2017.XIII.071.

[7] M. de Sa, P. Kotaru, K. Sreenath, "Point Cloud-Based Control Barrier Function Regression for Safe and Efficient Vision-Based Control," in *Proc. IEEE ICRA*, 2024. DOI: 10.1109/ICRA57147.2024.10610647.

### 신경망 음함수 계열

[8] J. Ortiz, A. Clegg, J. Dong, E. Sucar, D. Novotny, M. Zollhoefer, M. Mukadam, "iSDF: Real-Time Neural Signed Distance Fields for Robot Perception," in *Proc. RSS*, 2022. arXiv: 2204.02296.

[9] T. Chen, M. Gandhi, C. Culbertson, K. Iancu, M. Schwager, "CATNIPS: Collision Avoidance Through Neural Implicit Probabilistic Scenes," *IEEE Trans. Robotics*, vol. 40, pp. 2712–2733, 2024. DOI: 10.1109/TRO.2024.3387428.

[10] M. Pantic, C. Cadena, R. Siegwart, L. Ott, "Sampling-Free Obstacle Gradients and Reactive Planning in Neural Radiance Fields," arXiv: 2205.01389, 2022.

[11] Y. Zhou, Y. Shi, H. Mao, Z. Wang, Y. Meng, K. Wang, J. Lei, "Spacecraft Safe Robust Control Using Implicit Neural Representation for Geometrically Complex Targets in Proximity Operations," *IEEE Trans. Aerospace and Electronic Systems*, 2025. arXiv: 2507.13672.

### 우주 도메인 — 장애물 표현 및 GNC

[12] A. Badawy, C. R. McInnes, "On-Orbit Assembly Using Superquadric Potential Fields," *Journal of Guidance, Control, and Dynamics*, vol. 31, no. 1, pp. 30–43, 2008. DOI: 10.2514/1.28865.

[13] X. Chen, Z. Bai, H. Wang, X. Chen, H. Zhao, G. Sheng, "Obstacle Avoidance for Non-Cooperative Target Spacecraft with Gaussian Mixture Model," *Advances in Space Research*, vol. 68, no. 10, pp. 4217–4233, 2021. DOI: 10.1016/j.asr.2021.08.009.

[14] C. Cao, J. Yue et al., "Obstacle Avoidance for Spacecraft in Close Proximity Operations Using Gaussian Mixture Model and Fixed-Time Control," *Applied Sciences*, vol. 10, no. 17, p. 5986, 2020. DOI: 10.3390/app10175986.

[15] M. Leomanni, G. Bianchini, A. Garulli, A. Giannitrapani, R. Quartullo, "Explicit Model Predictive Control for Low-Thrust Spacecraft Proximity Operations," *Journal of Guidance, Control, and Dynamics*, vol. 45, no. 5, 2022. DOI: 10.2514/1.G006321.

[16] T. H. Park, S. D'Amico, "Rapid Abstraction of Spacecraft 3D Structure from Single 2D Image," in *AIAA SciTech Forum*, 2024. DOI: 10.2514/6.2024-0963.

[17] A. Vela-Rincón, M. Zamaro, N. Ortiz-Gómez, "GNC Design for Proximity Operations with a Tumbling Target in Retrograde Orbit," *Acta Astronautica*, vol. 196, pp. 380–393, 2022. DOI: 10.1016/j.actaastro.2022.04.026.

[18] G. Zhao, S. Xu, Y. Bo, "LiDAR-Based Non-Cooperative Tumbling Spacecraft Pose Tracking by Fusing Depth Maps and Point Clouds," *Sensors*, vol. 18, no. 10, p. 3432, 2018. DOI: 10.3390/s18103432.

### 3D 점군 피팅 및 표면 재구성

[19] B. S. Morse, T. S. Yoo, P. Rheingans, D. T. Chen, K. R. Subramanian, "Interpolating Implicit Surfaces From Scattered Surface Data Using Compactly Supported Radial Basis Functions," in *Proc. SMI*, 2001. DOI: 10.1109/SMA.2001.923379.

[20] M. Kazhdan, M. Chuang, S. Rusinkiewicz, H. Hoppe, "Poisson Surface Reconstruction with Envelope Constraints," *Computer Graphics Forum (SGP)*, 2020. DOI: 10.1111/cgf.14068.

[21] B. Dai, R. Khorrambakht, P. Krishnamurthy, F. Khorrami, "Sailing Through Point Clouds: Safe Navigation Using Point Cloud-Based Control Barrier Functions," *IEEE Robotics and Automation Letters*, vol. 9, no. 9, pp. 7731–7738, 2024. DOI: 10.1109/LRA.2024.3431870.

[22] W. Wen, Y. Yang, Z. Wang, T. Ren, W. Xu, "Implicit Swept Volume SDF: Enabling Continuous Collision-Free Trajectory Generation for Arbitrary Shapes," *ACM Trans. Graphics (SIGGRAPH)*, vol. 43, no. 4, 2024. DOI: 10.1145/3658181.

---

*최종 업데이트: 2026년 5월*
