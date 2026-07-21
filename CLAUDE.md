# lidar-2603

RANSAC 기반 고차 곡면(Fischer basis, ridge 정규화) 피팅 + LiDAR 포인트클라우드 정합 연구 프로젝트 (논문 브랜치: `paper/2026taes`).

## 실행 환경 (중요)

**MATLAB은 반드시 R2025a로 실행할 것.** 시스템에 R2026a도 설치돼 있지만(`which matlab`이 R2026a를 잡음), 이 프로젝트는 R2025a 기준이며 R2026a에는 필요한 툴박스(Symbolic Math Toolbox 등) 라이선스가 없어 `homogeneFischerTerms`의 `syms`부터 실패한다.

- 실행 파일: `C:\Program Files\MATLAB\R2025a\bin\matlab.exe`
- 배치 실행 예: `& "C:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "run('experiments/sec3_1_basis/lcurve_demo.m')"`
- 필요 툴박스: Symbolic Math + Statistics.

## 참고 문서

작업 전에 `docs/` 폴더를 먼저 읽을 것. 논문 섹션(`sec3_1_basis`, `sec3_2_convergence`, `sec3_3_iss`)별 배경, 실험 설계, 목표가 정리되어 있음.

## 폴더 구조

- 루트 `.m` 파일들 — 핵심 솔버/유틸 (`PoliNavigationSolver3_*`, `homogeneFischerTerms`, `regressionFourthOrder` 등)
- `experiments/` — 논문 섹션별 실험 스크립트. 모든 스크립트 상단에서 `experiments/setup_paths.m` 실행
- `utils/` — 실험 공용 함수 (`chamferDistance`, `applyOcclusion`, `spectralCutting`)
- `docs/` — 논문 섹션 설명 md (git 추적, ignore 대상 아님)
- `DATA/` — 실측 LiDAR PLY 스캔 (git-ignore)
- `IMAGE/` — 루트에서 생성되는 fig/png 산출물 (git-ignore)
- `test_ridge/` — ridge 스윕 산출물 (git-ignore)
