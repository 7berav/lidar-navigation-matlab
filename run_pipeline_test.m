% 스펙트럴 분할 파이프라인 테스트 스크립트
% 판 2개 교차 씬 → build_graph → normalize(alpha) → embed → segment → evaluate
%
% 가이드 체크리스트 대응:
%   진단: eigengap 프로파일 + Cheeger 상한 vs 실제 컷 conductance
%   컷:   MST-length(현행) vs MST-conductance(권장 최소 개선) 비교
rng(1);

%% 1) 합성 씬: 길쭉한 판 2개, 교차각 파라미터
aspect = 10;          % 판 aspect ratio (길이/두께 방향)
angDeg = 45;          % 교차각 [deg]
noise  = 0.02;        % 노이즈 sigma

Rz = @(th) [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
shapes = struct('type',{},'N',{},'size',{},'R',{},'t',{});
shapes(1) = struct('type','box','N',4000,'size',[aspect 2 0.3], ...
                   'R',eye(3),          't',[0 0 0]);
shapes(2) = struct('type','box','N',4000,'size',[aspect 2 0.3], ...
                   'R',Rz(deg2rad(angDeg)), 't',[0 0 1.2]);

[P, labels_gt] = synthPointCloud(shapes, noise);
fprintf('점군 %d개 생성 (aspect=%g, angle=%g deg, noise=%g)\n', ...
        size(P,1), aspect, angDeg, noise);

%% 2) 그래프 + 정규화 + 임베딩
gOpts = struct('k',30, 'gamma',1.0, 'minCompSize',150);
[W, gInfo] = buildGraph(P, gOpts);
P2  = P(gInfo.idxKeep, :);
gt2 = labels_gt(gInfo.idxKeep);

alpha = 1;                                % {0, 0.5, 1} 스윕 축
[L, DinvSqrt] = normalizeGraph(W, alpha);

kEig = 45;
[U, lambda] = embedSpectral(L, kEig);
Y = DinvSqrt * U;                         % random-walk 좌표

% --- 진단 1: eigengap 프로파일 ---
figure(10); clf;
subplot(1,2,1);
plot(1:kEig, lambda, 'o-', 'LineWidth', 1.2);
xlabel('index k'); ylabel('\lambda_k'); title('L spectrum'); grid on;
subplot(1,2,2);
plot(1:kEig-1, diff(lambda), 's-', 'LineWidth', 1.2);
xlabel('k'); ylabel('\lambda_{k+1}-\lambda_k'); title('eigengap'); grid on;

% --- 진단 2: Cheeger 상한 ---
lam2 = lambda(2);
fprintf('lambda_2 = %.3e  →  Cheeger 상한 sqrt(2*l2) = %.3e\n', ...
        lam2, sqrt(2*lam2));

%% 3) 분할: length(현행) vs conductance(개선) 비교
% eigengap으로 K 자동 추정: K* = argmax_k (lambda_{k+1}-lambda_k)
[Kgap, gapsAll] = selectEigenGap(lambda, struct('kMin',2, 'kMax',30));
fprintf('eigengap 추정 K* = %d (gap = %.3e)\n', Kgap, gapsAll(Kgap));

nSkip = 1;                                % 상수 고유벡터 건너뛰기
nDim  = max(Kgap - nSkip, 2);             % u_2..u_{K*} 사용 (수동 덮어쓰기 가능)
Yuse  = Y(:, nSkip+1 : nSkip+nDim);
maxCuts = max(Kgap - 1, 1);               % K* 군집이면 컷은 최대 K*-1개

sOptsL = struct('ky',40, 'cutMethod','length',      'qThr',0.96, ...
                'minSize',60, 'maxCuts',maxCuts);
sOptsC = struct('ky',40, 'cutMethod','conductance', 'qThr',0.90, ...
                'phiMax',0.05, 'minSize',60, 'maxCuts',maxCuts);

[labL, dgL] = segmentEmbedding(Yuse, W, sOptsL);
[labC, dgC] = segmentEmbedding(Yuse, W, sOptsC);

mL = evalSegmentation(labL, gt2);
mC = evalSegmentation(labC, gt2);
fprintf('\n%-14s  ARI=%.3f  K=%d (gt %d)  coverage=%.2f  cuts=%d\n', ...
    'length',      mL.ARI, mL.K_pred, mL.K_gt, mL.coverage, size(dgL.cutEdges,1));
fprintf('%-14s  ARI=%.3f  K=%d (gt %d)  coverage=%.2f  cuts=%d\n', ...
    'conductance', mC.ARI, mC.K_pred, mC.K_gt, mC.coverage, size(dgC.cutEdges,1));

% --- 참조 베이스라인: Fiedler 부호컷 (K=2) ---
% lambda_2가 작으면 2번째 고유벡터 부호만으로 2분할이 거의 완벽해야 함
fied = Y(:,2);
labF = uint32(1 + double(fied >= 0));
mF = evalSegmentation(labF, gt2);
fprintf('%-14s  ARI=%.3f  K=%d (gt %d)  coverage=%.2f\n', ...
    'fiedler-sign', mF.ARI, mF.K_pred, mF.K_gt, mF.coverage);

% conductance 채택된 컷들의 phi vs Cheeger 상한 (진단용, 실패해도 무시)
try
    phiC = dgC.cutEdges(:,4);
    fprintf('채택 컷 phi: min=%.3e  max=%.3e  (Cheeger 상한 %.3e)\n', ...
        min(phiC), max(phiC), sqrt(2*lam2));
catch ME
    fprintf('(phi 진단 출력 건너뜀: %s)\n', ME.message);
end

%% 4) 시각화: 정답 / Fiedler / length / conductance
titles = {'Ground truth', ...
          sprintf('Fiedler-sign (ARI %.2f)', mF.ARI), ...
          sprintf('MST-length (ARI %.2f)', mL.ARI), ...
          sprintf('MST-conductance (ARI %.2f)', mC.ARI)};
labSet = {gt2, labF, labL, labC};

figure(11); clf;
for s = 1:4
    subplot(2,2,s); hold on;
    lab = labSet{s};
    idx0 = (lab == 0);
    if any(idx0)
        scatter3(P2(idx0,1), P2(idx0,2), P2(idx0,3), 1, [0.6 0.5 0.5], ...
                 'filled', 'MarkerFaceAlpha', 0.4);
    end
    Ks = double(max(lab));
    cmap = lines(max(Ks,1));
    for kk = 1:Ks
        m = (lab == kk);
        if any(m)
            scatter3(P2(m,1), P2(m,2), P2(m,3), 3, cmap(kk,:), 'filled', ...
                     'MarkerFaceAlpha', 0.7);
        end
    end
    hold off; axis equal; grid on; view([1 1 1]);
    title(titles{s}); xlabel('X'); ylabel('Y'); zlabel('Z');
end
