# 3D 점군 피팅·표현 방식별 주요 논문 조사

**조사 기준:** 2D+높이 기반(BEV, elevation map)이 아닌 진정한 3D 점군(x, y, z 모두 독립적으로 처리)을 입력으로 사용하는 논문 한정.
**위성/우주 도메인 논문은 🛰 별도 표기.**

---

## A. 음함수 다항식 (Implicit Polynomial) 계열

### A1. 3L 알고리즘 — IP 피팅의 기준선

| 항목 | 내용 |
|---|---|
| **저자** | M. M. Blane, T. Lei, H. Çivi, D. B. Cooper |
| **제목** | The 3L Algorithm for Fitting Implicit Polynomial Curves and Surfaces to Data |
| **게재** | IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 22, no. 3, pp. 298–313 |
| **연도** | 2000 |
| **DOI** | 10.1109/34.841760 |

**3D 데이터 특징:** Princeton Shape Benchmark 3D 메시 샘플링 점군 (의자·화병·토끼 등 정적 객체).

**개선점:** 경계면 내·외부에 두 개 오프셋 수준집합을 구성하여 IP 피팅을 선형 시스템으로 정식화. 대수적 방법의 실용적 기준선 확립. 빠른 계산, 쉬운 구현이 장점이나 최적화 과정에서 실제 유클리드 거리 의미 부재로 정확도 한계 존재.

---

### A2. Ridge Regression에 의한 IP 안정화

| 항목 | 내용 |
|---|---|
| **저자** | T. Tasdizen, J.-P. Tarel, D. B. Cooper |
| **제목** | Improving the Stability of Algebraic Curves for Applications |
| **게재** | IEEE Transactions on Image Processing, vol. 9, no. 3, pp. 405–416 |
| **연도** | 2000 |
| **DOI** | 10.1109/83.826778 |

**3D 데이터 특징:** 3D 합성·실측 메시 점군.

**개선점:** 유클리드 불변 3D Ridge 행렬을 명시적 공식으로 유도하여 선형 피팅의 전역 불안정성 해소. 고차 단항식 기저의 조건수 문제를 사후(post-hoc) 페널티 항으로 완화. 단순 단위 행렬 Ridge 대비 좌표계 방향에 무관한 일관된 정규화 보장.

---

### A3. 제약된 IP — 영집합 형태 보장

| 항목 | 내용 |
|---|---|
| **저자** | D. Keren, C. Gotsman |
| **제목** | Fitting Curves and Surfaces With Constrained Implicit Polynomials |
| **게재** | IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 21, no. 1, pp. 31–41 |
| **연도** | 1999 |
| **DOI** | 10.1109/34.745731 |

**3D 데이터 특징:** 3D 물체 메시 샘플 점군.

**개선점:** 영집합이 star-shaped 또는 유계(bounded)인 곡면이 되도록 형태 위상학적 제약 부과. 다항식이 데이터 분포에 무관하게 닫힌 곡면을 생성하도록 보장. 위상 불일치(phantom sheets) 방지.

---

### A4. 음함수 B-spline (IBS) — 3L의 국소 기저 확장

| 항목 | 내용 |
|---|---|
| **저자** | M. Rouhani, A. D. Sappa |
| **제목** | Implicit B-Spline Fitting Using the 3L Algorithm |
| **게재** | Proceedings of IEEE International Conference on Computer Vision and Pattern Recognition (CVPR) |
| **연도** | 2010 |
| **비고** | MATLAB Central 구현 코드 공개 (fileexchange/44653) |

**3D 데이터 특징:** 3D LiDAR 포함 다양한 점군 (Stanford Bunny, 실측 물체 스캔).

**개선점:** 3L 선형 시스템 구조를 유지하면서 전역 단항식 기저를 격자 기반 B-spline 국소 기저로 교체. 설계 행렬이 구조적으로 희소해져 조건수 개선. 격자 해상도로 표현 유연성을 차수와 독립적으로 조절 가능. tension 정규화로 전체 곡면 형태 제어.

---

### A5. 적응형 안정 IP 피팅 — 기하적 방법 개선

| 항목 | 내용 |
|---|---|
| **저자** | J. Zheng, J. Takamatsu, K. Ikeuchi |
| **제목** | An Adaptive and Stable Method for Fitting Implicit Polynomial Curves and Surfaces |
| **게재** | IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 32, no. 3, pp. 561–568 |
| **연도** | 2010 |
| **DOI** | 10.1109/TPAMI.2009.189 |

**3D 데이터 특징:** 노이즈 포함 실측 3D 점군.

**개선점:** 직교 거리(orthogonal distance) 최소화 기반 기하적 방법. Levenberg-Marquardt 솔버 수렴 안정성 개선. 대수적 방법 대비 실제 유클리드 거리 의미 보존으로 정확도 향상. 단, 계산 비용이 높아 실시간 적용 제한.

---

### A6. SOS 다항식 + SDP 거리 계산 ★ 본 연구 직접 선조

| 항목 | 내용 |
|---|---|
| **저자** | A. A. Ahmadi, G. Hall, A. Makadia, V. Sindhwani |
| **제목** | Geometry of 3D Environments and Sum of Squares Polynomials |
| **게재** | Proceedings of Robotics: Science and Systems XIII (RSS) |
| **연도** | 2017, Cambridge, MA |
| **DOI** | 10.15607/RSS.2017.XIII.071 |
| **arXiv** | 1611.07369 |

**3D 데이터 특징:** Princeton Shape Benchmark 3D 메시 샘플링 점군 (정적 객체, 오프라인 처리).

**개선점:** SOS(Sum-of-Squares) 부수준집합을 이용해 두 객체 간 유클리드 거리, growth distance, penetration depth를 SDP(Semidefinite Program)로 계산. 경로계획기에 삽입 가능한 폐쇄형 거리 oracle 제공. 차수 2~8까지 적용.

**한계 (본 연구와의 차이):** 정적 메시 가정 (텀블링 타겟 불가), 이상치 처리 없음, 다물체 분리 없음, 우주 도메인 미적용.

---

### A7. 점군 기반 이차 CBF 회귀 ★ 본 연구 직접 비교군

| 항목 | 내용 |
|---|---|
| **저자** | M. de Sa, P. Kotaru, K. Sreenath |
| **제목** | Point Cloud-Based Control Barrier Function Regression for Safe and Efficient Vision-Based Control |
| **게재** | IEEE International Conference on Robotics and Automation (ICRA) |
| **연도** | 2024, Yokohama, Japan |
| **DOI** | 10.1109/ICRA57147.2024.10610647 |

**3D 데이터 특징:** Intel RealSense D435i 깊이 카메라 실시간 3D 점군 (매 제어 주기 갱신, 수천 점 규모).

**개선점:** 매 주기 국소 이차 다항식 h(q) = qᵀAq + bᵀq + c를 OLS로 피팅 → CBF-QP 직접 투입. SLAM 없이 점군에서 실시간 안전 보장. TurtleBot3(지상), Crazyflie(공중) 실험 검증.

**한계:** 차수 2 이차식으로 복잡 형상 표현 불가, 이상치 명시 처리 없음, 우주 도메인 미적용.

---

### A8. 다항식 분리 초곡면 (최신 프리프린트)

| 항목 | 내용 |
|---|---|
| **저자** | Y. Li et al. |
| **제목** | Online Trajectory Optimization for Arbitrary-Shaped Mobile Robots via Polynomial Separating Hypersurfaces |
| **게재** | arXiv 프리프린트 |
| **연도** | 2026 |
| **arXiv** | 2601.09231 |
| **비고** | 동료심사 전 프리프린트 |

**3D 데이터 특징:** 3D 이동로봇 환경 점군.

**개선점:** 로봇과 장애물 사이 다항식 분리 초곡면을 NLP 공동 최적화로 추정. 임의 비볼록 형상 로봇 적용 가능. 온라인 궤적 최적화와 직접 결합.

---

## B. 신경망 SDF / NeRF 계열

### B1. iSDF — 실시간 신경망 SDF 재구성

| 항목 | 내용 |
|---|---|
| **저자** | J. Ortiz, A. Clegg, J. Dong, E. Sucar, D. Novotny, M. Zollhoefer, M. Mukadam |
| **제목** | iSDF: Real-Time Neural Signed Distance Fields for Robot Perception |
| **게재** | Robotics: Science and Systems (RSS) |
| **연도** | 2022 |
| **arXiv** | 2204.02296 |
| **코드** | github.com/facebookresearch/iSDF |

**3D 데이터 특징:** Depth 카메라 실시간 3D 점군 스트림 (posed depth images).

**개선점:** MLP가 3D 좌표 → SDF 값을 온라인 연속 학습. 자기지도(self-supervised) 손실로 SDF 근사. Voxel Grid 대비 적응적 해상도, 컴팩트 표현, 부분 관측 영역 플로스인(fill-in). 충돌 비용·그래디언트를 직접 플래너에 제공.

---

### B2. CATNIPS — NeRF 기반 확률적 충돌 회피

| 항목 | 내용 |
|---|---|
| **저자** | T. Chen, M. Gandhi, C. Culbertson, K. Iancu, M. Schwager |
| **제목** | CATNIPS: Collision Avoidance Through Neural Implicit Probabilistic Scenes |
| **게재** | IEEE Transactions on Robotics, vol. 40, pp. 2712–2733 |
| **연도** | 2024 |
| **DOI** | 10.1109/TRO.2024.3387428 |
| **arXiv** | 2302.12931 |

**3D 데이터 특징:** RGB 카메라 입력 → NeRF 암묵적 3D 씬 (포즈 없는 단안 카메라만으로 훈련 가능).

**개선점:** NeRF를 Poisson Point Process(PPP)로 수학적 변환 → 충돌 확률 정량화. PURR(Probabilistically Unsafe Robot Region) 표현 도입. 궤적 계획에 확률적 안전 보장. NASA ULI 지원.

---

### B3. 신경망 SDF + CBF — 우주 도메인 SOTA 🛰 ★ 본 연구 직접 경쟁자

| 항목 | 내용 |
|---|---|
| **저자** | Y. Zhou, Y. Shi, H. Mao, Z. Wang, Y. Meng, X. Wang, L. Lei |
| **제목** | Spacecraft Safe Robust Control Using Implicit Neural Representation for Geometrically Complex Targets in Proximity Operations |
| **게재** | IEEE Transactions on Aerospace and Electronic Systems |
| **연도** | 2025 |
| **arXiv** | 2507.13672 |
| **비고** | 게재 예정, DOI 확인 필요 |

**3D 데이터 특징:** 비협조 우주 타겟 3D LiDAR 점군에서 직접 SDF 학습.

**개선점:** Enhanced implicit geometric regularization으로 보수적(over-approximation) 경계 학습. 2층 계층적 안전 강건 제어 (CBF 기반). Local minimum 완화를 위한 순환 부등식 도입. 비협조 우주 타겟 복잡 형상 대응.

**한계 (본 연구와의 차이):** MLP 파라미터 수만~수십만 개 (vs 본 연구 35~84개), 폐쇄형 그래디언트 없어 SDP 직접 삽입 불가, 단일 객체 가정.

---

### B4. 연속 암묵적 SDF 기반 임의 형상 궤적 최적화

| 항목 | 내용 |
|---|---|
| **저자** | S. Chen, P. Lu |
| **제목** | Continuous Implicit SDF Based Any-shape Robot Trajectory Optimization |
| **게재** | arXiv 프리프린트 |
| **연도** | 2023 |
| **arXiv** | 2303.01330 |

**3D 데이터 특징:** 임의 형상 로봇 주변 3D 점군.

**개선점:** 임의 형상 로봇(비구형)의 SDF를 연속 암묵적 표현으로 모델링하여 궤적 최적화에 직접 통합. 로봇 자체 기하 형상의 3D 표현 포함.

---

### B5. RMMI — 신경 SDF 기반 반응 조작 제어

| 항목 | 내용 |
|---|---|
| **저자** | (RMMI 저자진) |
| **제목** | RMMI: Enhanced Obstacle Avoidance for Reactive Mobile Manipulation using an Implicit Neural Map |
| **게재** | arXiv 프리프린트 |
| **연도** | 2024 |
| **arXiv** | 2408.16206 |

**3D 데이터 특징:** 3D 점군 → 신경 SDF 실시간 구성.

**개선점:** 신경 SDF 쿼리로 충돌 거리·방향 얻어 반응 제어기 부등식 제약 직접 삽입. 고정 쿼리 시간으로 일정한 제어 주기 보장. 로봇 표면 샘플 포인트 기반 근접성 평가.

---

## C. 기하 프리미티브 피팅 계열 (초이차·타원체·GMM)

### C1. Multi-scale 초이차 피팅 — 3D 점군 자세 추정

| 항목 | 내용 |
|---|---|
| **저자** | K. Duncan et al. |
| **제목** | Multi-scale Superquadric Fitting for Efficient Shape and Pose Recovery of Unknown Objects |
| **게재** | IEEE International Conference on Robotics and Automation (ICRA) |
| **연도** | 2013 |
| **URL** | cse.usf.edu/~sarkar/PDFs/DuncanICRA2013.pdf |

**3D 데이터 특징:** RGB-D (Kinect) 3D 점군.

**개선점:** 다중 스케일 복셀화로 초이차 파라미터(형상, 크기, 자세) 빠른 추정. 볼륨 기반 부품 표현으로 grasping 자세 추출. 삼축 대칭 가정으로 가정용 물체 근사.

---

### C2. 확률적 초이차 복원 — 이상치 강건 CVPR 2022

| 항목 | 내용 |
|---|---|
| **저자** | L. Liu et al. |
| **제목** | Robust and Accurate Superquadric Recovery: A Probabilistic Approach |
| **게재** | IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR) |
| **연도** | 2022 |
| **URL** | openaccess.thecvf.com/content/CVPR2022/papers/Liu_Robust_and_Accurate_Superquadric_Recovery_A_Probabilistic_Approach_CVPR_2022_paper.pdf |

**3D 데이터 특징:** 3D 점군 (ShapeNet, ScanNet 등 다양한 실측·합성 데이터).

**개선점:** GMM + EM 프레임워크로 초이차 피팅을 MLE 문제로 정식화. 이상치 오염 60%까지 강건 (하이퍼파라미터 튜닝 불필요). 축 비율 제한 없음. 계층적 다중 초이차 복원 가능.

---

### C3. 2D 이미지 → 3D 초이차 조립 — 우주 도메인 🛰

| 항목 | 내용 |
|---|---|
| **저자** | T. H. Park, S. D'Amico |
| **제목** | Rapid Abstraction of Spacecraft 3D Structure from Single 2D Image |
| **게재** | AIAA SciTech Forum |
| **연도** | 2024, Orlando, FL |
| **DOI** | 10.2514/6.2024-0963 |

**3D 데이터 특징:** 단일 2D 이미지에서 CNN으로 3D 초이차 구조 예측 (3D 점군 직접 피팅 아님).

**개선점:** CNN이 초이차 어셈블리의 형태 파라미터, 크기, 자세를 단일 이미지에서 직접 예측. 우주 타겟의 컴팩트 3D 구조 표현. NeRF 대비 낮은 데이터셋 요구량.

**한계:** 3D 점군 직접 피팅 아님, CAD/학습 데이터 의존.

---

### C4. GMM-APF — 복잡 형상 비협조 타겟 우주 도메인 🛰

| 항목 | 내용 |
|---|---|
| **저자** | X. Chen, Z. Bai, H. Wang, X. Chen, H. Zhao, G. Sheng |
| **제목** | Obstacle Avoidance for Non-Cooperative Target Spacecraft with Gaussian Mixture Model |
| **게재** | Advances in Space Research, vol. 68, no. 10, pp. 4217–4233 |
| **연도** | 2021 |
| **DOI** | 10.1016/j.asr.2021.08.009 |

**3D 데이터 특징:** 비협조 타겟 3D CAD 모델 샘플링 (점군 직접 측정 아님).

**개선점:** GMM 형태로 척력 포텐셜(APF) 구성. 회전 타겟 대응 최초 GMM 적용. Equal-collision-probability-curve로 위치별 충돌 확률 정량화.

---

### C5. GMM-APF + 고정시간 제어 우주 도메인 🛰

| 항목 | 내용 |
|---|---|
| **저자** | C. Cao, J. Yue et al. |
| **제목** | Obstacle Avoidance for Spacecraft in Close Proximity Operations Using Gaussian Mixture Model and Fixed-Time Control |
| **게재** | Applied Sciences, vol. 10, no. 17, p. 5986 |
| **연도** | 2020 |
| **DOI** | 10.3390/app10175986 |

**3D 데이터 특징:** 텀블링 타겟 3D 형상 샘플.

**개선점:** GMM-APF에 고정시간(fixed-time) 수렴 제어 결합. 임무 시간 제약 만족. 텀블링 비협조 타겟 우주 근접 운용 적용.

---

### C6. 초이차 + Voronoi 방향 포함 경로계획

| 항목 | 내용 |
|---|---|
| **저자** | L. Yang, G. Iyer, B. Lou, S. H. Turlapati, C. Lv, D. Campolo |
| **제목** | Path Planning in Complex Environments with Superquadrics and Voronoi-Based Orientation |
| **게재** | arXiv 프리프린트 |
| **연도** | 2024 |
| **arXiv** | 2411.05279 |

**3D 데이터 특징:** 3D 환경 점군.

**개선점:** 초이차 표현과 Voronoi 다이어그램 결합. 방향(orientation) 포함 복잡 3D 환경 경로계획. 다물체 처리 가능.

---

### C7. 강건 타원체 피팅 — EM 기반

| 항목 | 내용 |
|---|---|
| **저자** | 저자진 (arXiv) |
| **제목** | Robust Ellipsoid-specific Fitting via Expectation Maximization |
| **게재** | arXiv 프리프린트 |
| **연도** | 2021 |
| **arXiv** | 2110.13337 |

**3D 데이터 특징:** 3D 점군 (이상치 포함 합성·실측 데이터).

**개선점:** LS 원리 탈피. 단위 구 위 샘플 GMM으로 타원체 파라미터 MLE. 이상치 60%까지 강건. 축 비율 제한 없음. ε-알고리즘으로 EM 수렴 가속.

---

## D. Poisson · RBF · 고전 음함수 표면 재구성 계열

### D1. Compact Support RBF — 대규모 점군 실용화

| 항목 | 내용 |
|---|---|
| **저자** | B. S. Morse, T. S. Yoo, P. Rheingans, D. T. Chen, K. R. Subramanian |
| **제목** | Interpolating Implicit Surfaces From Scattered Surface Data Using Compactly Supported Radial Basis Functions |
| **게재** | Proceedings of Shape Modeling International (SMI) |
| **연도** | 2001 |
| **DOI** | 10.1109/SMA.2001.923379 |

**3D 데이터 특징:** 분산 3D 표면 점군 (Stanford Bunny 최대 35,947점).

**개선점:** Thin-plate spline(전역 지지, O(n²)) → Wendland Compact Support RBF로 교체. 행렬 희소화로 복잡도 O(n log n), 메모리 O(n)으로 감소. k-d tree 기반 근방 탐색. 등위면 추출 및 미분 기하 분석 가능.

---

### D2. Screened Poisson + 봉투 제약

| 항목 | 내용 |
|---|---|
| **저자** | M. Kazhdan, M. Chuang, S. Rusinkiewicz, H. Hoppe |
| **제목** | Poisson Surface Reconstruction with Envelope Constraints |
| **게재** | Computer Graphics Forum (SGP) |
| **연도** | 2020 |
| **DOI** | 10.1111/cgf.14068 |

**3D 데이터 특징:** RGB-D 3D 스캔 점군 (결손 데이터 포함).

**개선점:** 시각적 헐(visual hull) / 깊이 헐(depth hull)로 Dirichlet 경계 조건 부과. 결손 영역 재구성 품질 향상. 표준 Screened Poisson 대비 비워터타이트(non-watertight) 경계 처리 개선.

---

### D3. p-Poisson + Curl-free — 법선 없는 재구성

| 항목 | 내용 |
|---|---|
| **저자** | M. Kang et al. |
| **제목** | p-Poisson Surface Reconstruction in Curl-free Flow from Point Clouds |
| **게재** | Advances in Neural Information Processing Systems (NeurIPS) |
| **연도** | 2023 |
| **URL** | proceedings.neurips.cc/paper_files/paper/2023/file/bd18189308a4c45c7d71ca83acf3deaa-Paper-Conference.pdf |

**3D 데이터 특징:** 법선 없는 원시(raw) 3D 점군.

**개선점:** 비선형 p-Poisson 방정식 + curl-free 조건 도입. Eikonal 방정식 기반 방법의 비유일해 문제 해소. 표면 법선 없이도 세밀한 형상 복원. Deep learning 기반 INR과 결합 가능.

---

### D4. 공간 분할 국소 IP 피팅

| 항목 | 내용 |
|---|---|
| **저자** | Y. Li, H. Borouchaki, H. Miao, J. Zhang |
| **제목** | Implicit Function-based 3D Reconstruction for Point Cloud Data |
| **게재** | Annals of Mathematics and Physics, vol. 8, no. 5, pp. 202–208 |
| **연도** | 2025 |
| **DOI** | 10.17352/amp.000164 |

**3D 데이터 특징:** 3D 점군 (공간 분할 Octree 기반).

**개선점:** Octree 기반 공간 분할 후 국소 음함수 피팅. 법선 추정·방향 일관성을 위한 하이브리드 방법 제안. 대규모 복잡 형상 처리에 적용 가능.

---

## E. 🛰 우주 도메인 — 3D 점군 특화

### E1. 텀블링 타겟 LiDAR 자세 추적 🛰

| 항목 | 내용 |
|---|---|
| **저자** | G. Zhao, S. Xu, Y. Bo |
| **제목** | LiDAR-Based Non-Cooperative Tumbling Spacecraft Pose Tracking by Fusing Depth Maps and Point Clouds |
| **게재** | Sensors, vol. 18, no. 10, p. 3432 |
| **연도** | 2018 |
| **DOI** | 10.3390/s18103432 |

**3D 데이터 특징:** 비협조 타겟 3D LiDAR 점군 (텀블링 운동 수치 시뮬레이션).

**개선점:** ICP + 적응형 복셀 그리드 다운샘플링 → 실시간 6-DOF 상대 자세 추정. 큰 자세 변화 대응. 깊이 맵과 점군 융합으로 정확도 향상.

**용도:** 자세 추정 (표면 함수 피팅 아님).

---

### E2. 텀블링 타겟 GNC 설계 — 역행 궤도 🛰

| 항목 | 내용 |
|---|---|
| **저자** | A. Vela-Rincón, M. Zamaro, N. Ortiz-Gómez |
| **제목** | GNC Design for Proximity Operations with a Tumbling Target in Retrograde Orbit |
| **게재** | Acta Astronautica, vol. 196, pp. 380–393 |
| **연도** | 2022 |
| **DOI** | 10.1016/j.actaastro.2022.04.026 |

**3D 데이터 특징:** 텀블링 타겟 3D 메시 기반 형상 추정.

**개선점:** 역행 궤도 텀블링 타겟 대상 완전 GNC 루프 설계. 형상 추정 후 근접 운용 알고리즘 연결. 다목적 최적화 기반 접근 궤적 계획.

**용도:** GNC 설계 (장애물 표현은 메시 수준).

---

### E3. LiDAR 점군 기반 비협조 위성 자세 추정 — 대칭성 처리 🛰

| 항목 | 내용 |
|---|---|
| **저자** | L. Renaut, H. Frei, A. Nüchter |
| **제목** | CNN-based Pose Estimation of a Non-Cooperative Spacecraft with Symmetries from LiDAR Point Clouds |
| **게재** | IEEE Transactions on Aerospace and Electronic Systems, pp. 1–16 |
| **연도** | 2024 |
| **DOI** | 10.1109/TAES.2024 (확인 필요) |

**3D 데이터 특징:** 비협조 위성 3D LiDAR 점군 (위성 대칭성 명시 처리).

**개선점:** CNN으로 LiDAR 점군에서 직접 자세 추정. 위성 구조의 회전 대칭 모호성(symmetry ambiguity) 처리. 실측 LiDAR 스캔 검증.

**용도:** 자세 추정 (표면 함수 피팅 아님).

---

### E4. 다면체 KOZ 기반 저추력 MPC 우주 도메인 🛰

| 항목 | 내용 |
|---|---|
| **저자** | M. Leomanni, G. Bianchini, A. Garulli, A. Giannitrapani, R. Quartullo |
| **제목** | Explicit Model Predictive Control for Low-Thrust Spacecraft Proximity Operations |
| **게재** | Journal of Guidance, Control, and Dynamics, vol. 45, no. 5 |
| **연도** | 2022 |
| **DOI** | 10.2514/1.G006321 |
| **arXiv** | 2107.07254 |

**3D 데이터 특징:** 다면체 KOZ (Keep-Out Zone) 표현 — 3D 점군 직접 피팅 아님.

**개선점:** 가변 수평선 LP 기반 명시적 MPC. 다면체 근사로 비볼록 KOZ 처리. 저추력 위성 실시간 구현 가능.

---

## F. 자율주행·드론 도메인 — 3D 점군 기반 장애물 처리

### F1. 3D LiDAR 기반 동적 장애물 탐지·추적 (자율주행)

| 항목 | 내용 |
|---|---|
| **저자** | A. Saha, B. C. Dhara |
| **제목** | 3D LiDAR-Based Obstacle Detection and Tracking for Autonomous Navigation in Dynamic Environments |
| **게재** | International Journal of Intelligent Robotics and Applications, vol. 8, pp. 39–60 |
| **연도** | 2024 |
| **DOI** | 10.1007/s41315-023-00302-1 |

**3D 데이터 특징:** 3D LiDAR 점군 → u-depth / restricted v-depth 표현 (3D 처리 후 투영).

**개선점:** 3D LiDAR 점군에서 직접 장거리 장애물 추정. u-depth·v-depth 표현으로 탐지 효율 향상. 동적 환경 실시간 추적.

---

### F2. Implicit Swept Volume SDF (드론 궤적)

| 항목 | 내용 |
|---|---|
| **저자** | W. Wen, Y. Yang, Z. Wang, T. Ren, W. Xu |
| **제목** | Implicit Swept Volume SDF: Enabling Continuous Collision-Free Trajectory Generation for Arbitrary Shapes |
| **게재** | ACM Transactions on Graphics (SIGGRAPH), vol. 43, no. 4, Art. 110 |
| **연도** | 2024 |
| **DOI** | 10.1145/3658181 |

**3D 데이터 특징:** 임의 형상 로봇의 3D swept volume SDF (신경망 기반).

**개선점:** 임의 형상 로봇의 swept volume을 신경망 SDF로 연속 표현. 연속 충돌 없는 궤적 생성. 형상 복잡도에 무관한 계산 효율.

---

### F3. 점군 SOS 궤적 최적화 (지상 + 쿼드로터)

| 항목 | 내용 |
|---|---|
| **저자** | Y. Li, C. Zheng, K. Chen, Y. Xie, X. Tang, M. Y. Wang, J. Ma |
| **제목** | Collision-Free Trajectory Optimization in Cluttered Environments Using Sums-of-Squares Programming |
| **게재** | arXiv 프리프린트 |
| **연도** | 2024 |
| **arXiv** | 2404.05242 |

**3D 데이터 특징:** 시뮬레이션 + 실제 환경 점군 (IRIS 자유공간 추출 후 처리).

**개선점:** 로봇 자체 기하를 반대수집합으로 표현. 자유공간 내 contain 관계를 SOS-SDP로 검증. KKT 라그랑지안에서 그래디언트 추출 → NLP 변환. 4~7 자유도 로봇 적용.

---

### F4. 점군 기반 Sailing CBF (사족보행 로봇)

| 항목 | 내용 |
|---|---|
| **저자** | B. Dai, R. Khorrambakht, P. Krishnamurthy, F. Khorrami |
| **제목** | Sailing Through Point Clouds: Safe Navigation Using Point Cloud-Based Control Barrier Functions |
| **게재** | IEEE Robotics and Automation Letters, vol. 9, no. 9, pp. 7731–7738 |
| **연도** | 2024 |
| **DOI** | 10.1109/LRA.2024.3431870 |
| **arXiv** | 2403.18206 |

**3D 데이터 특징:** 3D 점군 직접 처리 (명시적 피팅 없음).

**개선점:** 점군에서 직접 스케일링 인수 추출 → CBF-QP. de Sa 2024의 후속 비교군. 피팅 단계 없이 점군 직접 활용.

---

## G. 비교 종합

### G1. 3D 점군 처리 방식별 핵심 트레이드오프

| 표현 방식 | 대표 논문 | 파라미터 수 | 해석적 그래디언트 | SDP 호환 | 이상치 강건 | 우주 도메인 |
|---|---|---|---|---|---|---|
| 음함수 다항식 (차수 4) | Blane 2000, Ahmadi 2017 | 35 | ✓ (폐쇄형) | ✓ | △ | ✗ (선행 연구 없음) |
| 음함수 다항식 (차수 6) | Tasdizen 2000 | 84 | ✓ (폐쇄형) | ✓ | △ | ✗ (선행 연구 없음) |
| 이차 IP (CBF 회귀) | de Sa 2024 | 10 | ✓ | ✓ | ✗ | ✗ |
| 초이차 | Duncan 2013, Park 2024 | 5–15 | ✓ | △ | △ | 🛰 일부 |
| GMM | Chen 2021, Cao 2020 | ~90 | ✓ | ✗ | △ | 🛰 일부 |
| 신경망 SDF | Ortiz 2022, Zhou 2025 | 수만~수십만 | △ (자동미분) | ✗ | △ | 🛰 Zhou 2025 |
| NeRF | CATNIPS 2024 | 수만~수십만 | △ | ✗ | △ | ✗ |
| **본 연구 (HERO)** | **정우진 2026** | **35~84** | **✓ (폐쇄형)** | **✓** | **✓ (RANSAC)** | **🛰 최초** |

### G2. 우주 도메인 3D 점군 피팅 갭 확인

우주 도메인 논문(Zhao 2018, Renaut 2024, Vela-Rincón 2022)은 모두 3D LiDAR 점군을 사용하나 **목적이 자세 추정**이며, 점군에서 해석적 장애물 표면 모델을 피팅하여 GNC 제약으로 연결하지 않는다. Zhou 2025(T-AES)가 유일하게 이 루프를 구성하나 신경망 SDF 기반이다.

**"3D 점군 → 비신경망 음함수 다항식 피팅 → GNC 제약"의 루프를 비협조 우주 도메인에서 완결한 동료심사 논문은 현재까지 존재하지 않는다.** 이것이 본 연구의 핵심 리서치 갭이다.

---

*최종 업데이트: 2026년 5월*
