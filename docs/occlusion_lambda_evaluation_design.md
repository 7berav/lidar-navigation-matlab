# Occlusion 환경 정규화(λ) 튜닝·평가 설계 정리

동차다항식 레벨셋 피팅(HERO 기반) 모델에 대해, 정규화 파라미터 λ를 튜닝하고 occlusion 강건성을 평가하기 위한 실험 설계 정리.

---

## 0. 전체 구조 (3단계)

| 단계 | 목적 | 핵심 산출물 |
|------|------|-------------|
| 1. Oracle 분석 | GT 기반 오차의 λ 의존성과 최적 basin 존재 확인 | 오차(SDF/CD/HD/IoU) vs λ의 U자 곡선 |
| 2. 실용적 선택 검증 | GT 없이 λ를 고를 수 있음을 입증 | L-curve corner(+GCV)가 oracle basin 안에 위치 |
| 3. 이식성·개선폭 검증 | 선택한 λ가 조건(occlusion/outlier) 넘어 유효함을 입증 | λ ∈ {0, λ*, ...}의 occlusion sweep |

핵심 메시지: "λ를 얼추 정했다"가 아니라 **"GT 없이 재현 가능한 λ 선택 절차를 제시하고 그 타당성을 검증했다"**로 서술.

---

## 1. 평가 지표

### 1.1 Chamfer Distance가 필수인가
- 필수는 아님. CD는 학습 기반 completion 커뮤니티(PCN, Yuan et al., 3DV 2018)의 관례적 표준일 뿐, 평가의 논리적 필요조건은 아님.
- 현재 쓰는 SDF 오차는 모델이 implicit level-set이라는 점과 직접 정합되므로 그 자체로 타당.
- 다만 **넣는 것을 권장**: (1) SDF 오차는 단방향이라 비관측 영역의 잉여 표면(spurious lobe)을 못 잡음 — 고차 다항식 과적합의 병리가 정확히 이 사각지대. 양방향 CD는 "빠진 형상"과 "생겨난 잉여 형상"을 모두 측정. (2) 타 연구와의 비교가능성(심사 지적 방어).
- CD는 평균 기반이라 outlier에 둔감·시각 품질과 괴리(Tatarchenko et al., CVPR 2019) → **F-score@τ 병행**이 최근 관례.

### 1.2 지표별 역할
| 지표 | 측정 대상 | 비고 |
|------|-----------|------|
| SDF 오차 | GT 표면점에서 다항식 값 오차 | 현재 사용 중. 단방향. |
| Chamfer Distance (양방향) | 점집합 간 최근접 거리의 평균(양방향 합) | 정확도·완전성이 한 수치로 뭉개짐 |
| **F-score@τ** | 허용오차 τ 내 올바르게 복원된 표면 비율 | P/R 분리. occlusion엔 recall이 핵심 |
| Hausdorff (HD) / max error | 최악 오차 | 충돌회피 안전여유와 직결. HD95 권장 |
| Volumetric IoU | 두 폐곡면 내부의 교집합/합집합 | **해석적 inside/outside라 사실상 공짜 — 모델 정합성 최상** |
| Normal Consistency | 대응점 법선 코사인 유사도 | 표면 기울기 오차 |

### 1.3 F-score@τ 정의
- Precision P(τ): 추정점 중 최근접 GT까지 거리가 τ 이내인 비율 → 잉여 형상 벌점.
- Recall R(τ): GT점 중 최근접 추정점까지 거리가 τ 이내인 비율 → 누락 형상 벌점.
- F(τ) = 2·P·R / (P+R), 범위 [0,1].
- CD 대비 장점: 정확도·완전성을 분리한 뒤 조화평균 → 한쪽만 나빠도 급락.
- **τ 관례**: 물체 특성크기의 1%. 본 세팅(1 m 물체, σ=0.01 m)이면 τ = 1 cm가 자연스럽고 노이즈·안전여유와 물리적으로 연결해 정당화 가능. 단일 τ 의존이 걱정되면 τ sweep 곡선으로 보고.

### 1.4 최소 권장 구성
양방향 CD + F-score@1%(P,R 분리) + Hausdorff(또는 max error) + volumetric IoU.
→ IoU는 폐곡면 레벨셋 특성을 가장 잘 살리므로 우선순위 높게.

---

## 2. GT 기준: CAD인가, occlusion 0% 추출 점군인가

**원칙: 완전한 GT 표면(CAD/해석적)에서 균일 재샘플링한 클린 점군을 기준으로.**

근거:
1. CD는 점 밀도 분포에 민감한데, 센서 시뮬 추출 점군은 시점 의존적으로 밀도 불균일.
2. 추출 점군엔 Gaussian 노이즈가 섞여 있어 기준으로 쓰면 센서 노이즈가 metric에 혼입 → 모델 오차와 분리 안 됨.
3. occlusion 강건성의 핵심은 비관측 영역 복원 품질 — 관측 점군 기준이면 그 영역이 평가에서 빠짐.
- (참고) PCN, Completion3D, MVP 모두 GT를 CAD 메시의 uniform/Poisson-disk 샘플(8k–16k점)로 정의.

**두 지표 병행 보고 권장**:
- CAD 기준 CD = reconstruction error(진짜 형상 복원도). ← 튜닝 목적함수에 해당.
- 관측 점군 기준 CD = data fidelity(주어진 관측 적합도).
- 가려진 영역만 따로 계산한 missing-region CD → occlusion 강건성 주장에 가장 직접적 근거.

---

## 3. Occlusion 절단 방식은 실제로 통용되는가

**네, 거리 순서 절단은 출판된 표준 프로토콜.**
- PoinTr(Yu et al., ICCV 2021) ShapeNet-55/34: 임의 시점 기준 먼 순서로 25/50/75% 제거 → easy/moderate/hard.
- PF-Net(Huang et al., CVPR 2020): 시점 기준 crop. seed point의 kNN 비율 제거(random patch removal)도 흔함.
- 3DLoMatch(CVPR 2021): overlap 10–30% 저중첩을 별도 벤치마크로.
- 현재 쓰는 z-quantile 절단은 이 방향성 절단 계열 → 관례상 문제없음.

**비현실적인데도 쓰는 이유**: 통제 가능성. occlusion을 단일 스칼라로 환원해야 sweep·Monte Carlo가 가능. ablation엔 오히려 정석.

**"보통 50%는 가려진다"에 대한 방어**:
- 단일 시점 스캔은 뒷면이 원리적으로 안 보임(self-occlusion) → 기본 가시율 이미 ≤50%.
- 그래서 completion 문헌은 가상 센서 렌더링으로 partial scan 생성(PCN: 8시점 depth 렌더 back-projection, hidden point removal(Katz et al., SIGGRAPH 2007), raycasting).
- **관례 이원화**: 통제 sweep엔 비율 절단, 실제성 주장엔 시점 기반 렌더링.

**제안**: quantile sweep은 Monte Carlo 통제 변수로 유지 + HPR/raycasting 단일시점 스캔(가시율 ~40–50%)을 별도 실제성 시나리오로 추가. RPOD 표적은 tumbling이므로 시간 경과에 따라 다중 시점 누적 → 유효 가시율 회복. "단일 프레임 ≤50%지만 회전 누적으로 유효 가시율 상승" 실험은 절단 비현실성 지적을 응용 적합성 논거로 전환.

---

## 4. L-curve의 위치와 역할

### 4.1 L-curve vs oracle 오차곡선은 다른 것
- **L-curve**(Hansen, 1992): 잔차 노름 ‖Ac−b‖ 대 페널티 노름 ‖Lc‖. **GT 없이** λ를 고르는 도구.
- **CD/HD/SDF vs λ**: GT를 아는 합성 실험에서만 그리는 **oracle 곡선**.
- 둘을 섞어 "L-curve로 λ 정했다"고 쓰면 "GT 있는데 왜 L-curve? / GT 없는 실전은?" 양쪽에 취약 → **반드시 분리 서술**.

### 4.2 서술 순서
1. **Oracle 분석**: 단일 조건에서 λ를 로그 sweep, GT 지표 vs λ의 U자. 보여야 할 것은 최적점보다 **평탄한 basin의 폭**(1–2 order면 "λ 선택에 둔감" = 강건성 주장).
2. **실용적 선택 검증**: L-curve corner가 oracle basin 안에 떨어짐 → "GT 없이 온라인 선택 가능" 근거. λ 선택을 절차가 아닌 **결과(result)**로 격상.
3. L-curve는 수렴 반례 존재(Vogel 1996, Hanke 1996) → **Morozov discrepancy(σ 기지) / GCV(Golub-Heath-Wahba 1979) 병행**해 두세 선택자가 같은 basin을 가리키면 논거 강화.

### 4.3 L-curve 잔차와 Chamfer/Hausdorff의 차이
| 축 | L-curve 잔차 | CD/HD |
|----|-------------|-------|
| 측정 공간 | 대수적(algebraic): 함수값 오차 \|f(x)−1\| | 기하학적(geometric): 유클리드 거리 |
| 관계 | 비례 안 함. d ≈ \|f−1\|/‖∇f‖ (Taubin, PAMI 1991). 고차일수록 괴리 심함 | — |
| 평가점 | **관측(가려진) 입력점만** → 비관측 영역·잉여 성분 못 봄. λ→0에서 단조 감소 | GT 기준 → U자, 과적합 탐지 |
| 역할 | 내부(internal), GT 불필요 → 실전 배치용 | 외부(external) 검증, 합성에서만 |

이 구분 명시가 "왜 두 종류 그래프가 다 필요한가"를 자연히 정당화.

---

## 5. 튜닝 시 outlier를 빼도 되는가

**빼는 것이 오히려 방법론적으로 깔끔.**
1. 역할 분리: λ는 ill-conditioning·과적합 담당, outlier 제거는 RANSAC 담당. 겹치면 효과 귀속 불가.
2. outlier 없으면 (노이즈 수준에서) 전 점이 사실상 inlier → 설계행렬 A·관측 b가 λ에 고정. **L-curve 이론은 고정 (A,b)에서 λ만 sweep을 전제**(Hansen 1992)하므로 outlier-free가 L-curve가 잘 정의되는 유일 세팅에 가까움.

**단, 검증 단계는 별도 필요**: λ가 RANSAC 루프 내 가설 적합에도 쓰이므로, clean에서 고른 λ*가 오염 하에서도 유효한지는 자명하지 않음. 기존 outlier sweep(2/5/10%) 인프라로 λ∈{0, λ*} end-to-end 비교. 실전 논리와 일치: 배치 시 L-curve는 RANSAC 통과 후 inlier 집합에서 계산되므로 outlier-free 튜닝은 "RANSAC 통과 후" 상태를 모사.

---

## 6. L-curve/CD/HD를 "어떤 결과치"로 그리는가

### 6.1 계수 c_λ는 어느 모델에서
**수렴까지 돌린 최종(재가중 포함) 모델.** 반복별 가설(minimal sample fit)은 부적합:
1. minimal sample 적합은 보간에 가까워 잔차가 구성상 ≈0 → 잔차 축 퇴화, corner 안 생김.
2. 가설은 분산 크고 재가중 전 중간 산물 → "고정 데이터에 대한 λ별 정칙화 해의 족"이 아님.
3. 튜닝 대상은 실제 배치할 추정기 = 수렴 최종 모델.

### 6.2 잔차를 어느 점 집합에서 평가
**최종 inlier 집합.** (outlier-free 조건에선 inlier 집합 ≈ 전체 점이므로 전체 점 사용 무방. 서술에 "outlier 없어 inlier=전체 점, 오염 검증 단계에선 inlier로 좁힘" 명기.)

### 6.3 CD/HD도 매 가설이 아니라 λ마다 최종 모델로
- CD/HD가 보려는 과적합은 **λ의 함수**이지 iteration의 함수 아님. 가로축이 λ인 곡선 → 각 λ마다 확정된 최종 결과 하나 필요.
- 가설 단위 CD/HD는 분산 폭발 + 다수는 버려질 후보 → 정규화·과적합 신호를 노이즈로 덮음. 과적합은 최종 추정기의 성질.
- U자(λ 클수록 과적합 억제로 하강 → 과하면 과소적합으로 상승)는 최종 모델 수준에서만 정의됨.

### 6.4 통합 프로토콜
outlier-free 조건에서 점군(노이즈 실현+가려짐) 하나 고정 → 로그 간격 λ마다 재가중 포함 전체 절차를 수렴까지 실행 → 최종 c_λ에서:
- (a) inlier 잔차 노름·페널티 노름 → L-curve 위 한 점
- (b) marching cubes로 표면점 샘플 → GT와 CD/HD 계산 → oracle 곡선 위 한 점
- λ sweep으로 두 곡선 동시 완성.
- RANSAC·노이즈 확률성 → Monte Carlo 시드별 반복, 곡선 겹쳐 그리거나 corner λ 분포(중앙값·사분위) 보고. corner λ 분포가 좁으면 재현성 근거.

---

## 7. RANSAC 반복 횟수 고정

- 이론 bound N ≥ log(1−p)/log(1−wˢ)는 "순수 inlier 샘플이 한 번은 뽑힐 확률 p"만 보장, best 선택 보장 아님 → 안전계수 곱해 크게 고정하는 것이 표준.
- **λ 실험에서 반복 횟수는 튜닝 대상이 아니라 통제 변수** → 충분히 크게 고정하고 명기. 그래야 λ 효과가 정규화 효과이지 수렴 미완 잔여 효과가 아님을 보장.
- 수렴 횟수 자체의 최적화는 별개 기여(석사논문 Future Work의 adaptive termination) → 효율 개선 절로 분리.
- 안전장치: "N을 2배로 늘려도 λ 곡선·corner 유의미하게 안 바뀜"(saturation) 확인을 각주/부록에.

### 실험 절차 주의
- L-curve 잔차는 최종 inlier 집합 기준(outlier 포함 잔차는 곡선 왜곡).
- 각 (λ, occlusion) 셀마다 Monte Carlo 반복 → 평균 ± 표준편차(또는 분위) 밴드.
- RANSAC 하이퍼파라미터(inlier threshold ε 등)는 λ sweep 동안 **고정**하고 명기. ε과 λ가 동시에 움직이면 효과 분리 불가.

---

## 8. 암시적 레벨셋에서 표면점 추출 (CD/HD 계산 실제)

추정 결과는 점집합이 아니라 f(x)=1 암시적 곡면 → 점 대 점 CD/HD를 쓰려면 곡면 위 점을 뽑아야 함. **전 공간을 잘라 탐색하지 않음.**

### 8.1 표준: Marching Cubes (Lorensen & Cline, SIGGRAPH 1987)
1. bounding box에 정규 격자(예: 128³), 각 격자점에서 f 값만 평가(해석적이라 매우 쌈).
2. 각 셀 8꼭짓점 부호 패턴으로 등위면이 지나는 셀에서만 삼각형 생성(전 공간 탐색 아님, 국소 처리).
3. 나온 메시에서 균일(면적가중/Poisson-disk) 점 샘플링 → 추정 점집합.
4. GT 클린 점집합과 CD/HD 계산.
- 도구: `skimage.measure.marching_cubes` → 메시, `trimesh sample.sample_surface` → 표면점, `scipy.spatial.cKDTree` → 최근접 탐색(O(log n), 전수비교 아님).

### 8.2 정의
- CD(X,Y) = (1/|X|)Σ_{x∈X} min_{y∈Y}‖x−y‖² + (1/|Y|)Σ_{y∈Y} min_{x∈X}‖y−x‖²
- HD(X,Y) = max{ max_{x∈X} min_{y∈Y}‖x−y‖ , max_{y∈Y} min_{x∈X}‖y−x‖ }
- HD는 단일 최악점이 값을 지배 → **HD95(95백분위) 병행 권장**.

### 8.3 대안: Taubin 근사 거리
- GT 표면점마다 |f(x)−1|/‖∇f(x)‖로 근사 기하 거리 직접 계산 → marching cubes 없이 GT→추정면 방향(recall 성격, 가려진 영역 복원) 즉시 계산.
- 단 반대 방향(추정면→GT, precision, 잉여 표면)은 못 잡음 → 잉여 성분까지 보려면 marching cubes로 양방향.
- **병행**: recall 계열 빠른 확인엔 Taubin, precision 포함 정식 CD/HD엔 marching cubes.

---

## 9. 권장 핵심 Figure

1. **6차 모델 occlusion sweep, λ=0 vs λ=λ\* 오버레이**: 기존 6차가 occlusion 40%에서 SDF 0.35로 발산하던 것이 정규화로 억제되면 Future Work 가설(정규화가 고차 occlusion 민감성 완화)을 정면 입증. 4차는 악화만 안 되면 충분.
2. **오차 vs λ U자 + basin 표시** (oracle).
3. **L-curve + corner + (GCV/Morozov) 표시**, corner가 basin 안에 위치.
4. (여유 시) **(λ, occlusion) 오차 heatmap + occlusion별 oracle 최적 λ\*(q) 궤적**: λ\*(q) 거의 일정 → "한 λ로 전 구간 커버", 증가 추세 → "가려짐 심할수록 강한 정규화 필요". 어느 쪽이든 발견.
- λ 집합 표준: **{0, λ\*, λ\*/10, 10λ\*}**. **λ=0(무정규화) 필수 포함**(기존 기법 대비 baseline).

---

## 참고문헌 (본문 인용)
- Hansen (1992) — L-curve
- Golub, Heath & Wahba (1979) — GCV
- Vogel (1996), Hanke (1996) — L-curve 수렴 반례
- Taubin (PAMI 1991) — approximate distance
- Lorensen & Cline (SIGGRAPH 1987) — Marching Cubes
- Yuan et al. (3DV 2018) — PCN
- Tatarchenko et al. (CVPR 2019) — F-score, CD 비판
- Yu et al. (ICCV 2021) — PoinTr
- Huang et al. (CVPR 2020) — PF-Net
- Katz et al. (SIGGRAPH 2007) — Hidden Point Removal
- Fan et al. (CVPR 2017) — EMD
- Mescheder et al. (CVPR 2019) — Occupancy Networks, IoU
