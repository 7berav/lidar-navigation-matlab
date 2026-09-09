% 파라미터 스윕: {k, gamma, alpha, qThr, phiMax} 그리드를 합성 씬에서 ARI로 정량 비교
%
% 목적 (후보군 결정 근거 확보):
%   - alpha 스윕 {0, 0.5, 1} : Coifman–Lafon 정규화가 밀도 아티팩트를 얼마나
%     지우는지 확인. alpha=1(Laplace–Beltrami 극한)이 우세하면 보로노이/코탄
%     라플라시안(A2 후보)은 도입하지 않아도 된다는 근거가 됨.
%   - cutMethod length vs conductance : MST 절단 방식 비교 (B1 vs B2)
%   - eigengap 자동 K (selectEigenGap) 사용
%
% 출력: results table (ARI 내림차순), sweep_results.mat 저장

rng(1);

%% 씬: 교차 판 2장 + 원통 1개 (트러스/모듈 근사)
Rz = @(th) [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
shapes = struct('type',{},'N',{},'size',{},'R',{},'t',{});
shapes(1) = struct('type','box','N',4000,'size',[10 2 0.3], ...
                   'R',eye(3),              't',[0 0 0]);
shapes(2) = struct('type','box','N',4000,'size',[10 2 0.3], ...
                   'R',Rz(deg2rad(45)),     't',[0 0 1.2]);
shapes(3) = struct('type','cylinder','N',2500,'size',[0.8 0.8 4], ...
                   'R',eye(3),              't',[0 6 0]);
noise = 0.02;
[P, labels_gt] = synthPointCloud(shapes, noise);
fprintf('점군 %d개, 정답 K = %d\n', size(P,1), double(max(labels_gt)));

%% 스윕 그리드
grid_k      = [20 30];
grid_gamma  = [0.5 1.0];
grid_alpha  = [0 0.5 1];
grid_qThr   = [0.90 0.96];
grid_phiMax = [0.02 0.05 0.10];
methods     = {'length', 'conductance'};

kEig = 45;
rows = {};

for gk = grid_k
  for gg = grid_gamma
    % 그래프는 (k, gamma)마다 1회만 생성
    [W, gInfo] = buildGraph(P, struct('k',gk, 'gamma',gg, 'minCompSize',150));
    gt2 = labels_gt(gInfo.idxKeep);

    for ga = grid_alpha
      [L, DinvSqrt] = normalizeGraph(W, ga);
      [U, lambda]   = embedSpectral(L, kEig);
      Y = DinvSqrt * U;

      [Kgap, gapsAll] = selectEigenGap(lambda, struct('kMin',2, 'kMax',30));
      nSkip = 1;
      nDim  = max(Kgap - nSkip, 2);
      Yuse  = Y(:, nSkip+1 : nSkip+nDim);
      maxCuts = max(Kgap - 1, 1);

      for im = 1:numel(methods)
        method = methods{im};
        for gq = grid_qThr
          phiList = grid_phiMax;
          if strcmp(method, 'length'), phiList = NaN; end   % phiMax 무관
          for gp = phiList
            sOpts = struct('ky',40, 'cutMethod',method, 'qThr',gq, ...
                           'minSize',60, 'maxCuts',maxCuts);
            if ~isnan(gp), sOpts.phiMax = gp; end
            try
                [lab, ~] = segmentEmbedding(Yuse, W, sOpts);
                mtr = evalSegmentation(lab, gt2);
                rows(end+1,:) = {gk, gg, ga, string(method), gq, gp, ...
                                 Kgap, mtr.ARI, mtr.K_pred, mtr.K_gt, ...
                                 mtr.coverage}; %#ok<SAGROW>
                fprintf(['k=%2d g=%.1f a=%.1f %-11s q=%.2f phi=%5.2f | ' ...
                         'K*=%2d ARI=%.3f K=%d/%d cov=%.2f\n'], ...
                        gk, gg, ga, method, gq, gp, Kgap, mtr.ARI, ...
                        mtr.K_pred, mtr.K_gt, mtr.coverage);
            catch ME
                warning('스윕 실패 (k=%d g=%.1f a=%.1f %s q=%.2f phi=%.2f): %s', ...
                        gk, gg, ga, method, gq, gp, ME.message);
            end
          end
        end
      end
    end
  end
end

%% 결과 정리
results = cell2table(rows, 'VariableNames', ...
    {'k','gamma','alpha','cutMethod','qThr','phiMax', ...
     'Kgap','ARI','K_pred','K_gt','coverage'});
results = sortrows(results, 'ARI', 'descend');
disp(results(1:min(15,height(results)), :));

save('sweep_results.mat', 'results', 'shapes', 'noise');
fprintf('sweep_results.mat 저장 완료 (%d개 구성)\n', height(results));

%% 요약: alpha별 / cutMethod별 최고 ARI (A1 vs A2 판단 근거)
fprintf('\n[alpha별 최고 ARI]\n');
for ga = grid_alpha
    m = results.alpha == ga;
    fprintf('  alpha=%.1f : ARI=%.3f\n', ga, max(results.ARI(m)));
end
fprintf('[cutMethod별 최고 ARI]\n');
for im = 1:numel(methods)
    m = results.cutMethod == string(methods{im});
    fprintf('  %-11s : ARI=%.3f\n', methods{im}, max(results.ARI(m)));
end
