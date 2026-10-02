% ISS: 설정 간 차이가 난수 시드(thinning 순서)보다 큰지 확인
%
% thinning 순서가 바뀌면 점 집합이 달라지고 분할 점수도 달라진다. 한 시드에서의 점수 차이가
% 시드 간 표준편차보다 작으면 의미가 없으므로, 여러 시드의 평균과 표준편차로 비교한다.
% 시드마다 점 집합이 다르므로 정답(GLB 부품 라벨)도 매번 다시 옮긴다 (issPartLabels).
%
% 출력: out_cutting_ISS/seedcheck_<yyyyMMdd_HHmmss>/ summary.csv, results.mat

seeds = 1:6;
g0 = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',150);
P  = @(varargin) struct('ky',40, varargin{:});
%   이름, 법선 sigma(0=없음, -1=레거시 alpha 0), 모드 필터, 컷 옵션, 특수 규칙
cfg = {
  'C0_legacy',            -1,  0,  P('cutMethod','length','qThr',0,'minSize',60),                                                     'legacy'
  'noNormal_phi',          0,  60, P('cutMethod','phiorder','minSize',60,'candidates','all','phiMax',0.005),                          ''
  'C8a_old',               0.2, 0, P('cutMethod','phiorder','minSize',60,'candidates','length','lenQuantile',0.9),                    'cheeger'
  'lengthOnly_K16',        0.2, 60, P('cutMethod','lensize','minSize',60,'maxCuts',15),                                                ''
  'phiOrder',              0.2, 60, P('cutMethod','phiorder','minSize',60,'candidates','all','phiMax',0.005),                          ''
  'lengthOrder_phiGate',   0.2, 60, P('cutMethod','phiorder','minSize',60,'candidates','all','phiMax',0.005,'lenPower',8),             ''
  'C12_final',             0.3, 60, P('cutMethod','phiorder','minSize',100,'candidates','all','phiMax',0.005,'lenPower',1),            ''
};
nC = size(cfg,1); nS = numel(seeds);
ARI = nan(nC,nS); PUR = ARI; FRG = ARI; KK = ARI; COV = ARI;
for si = 1:nS
    P_use = prepareISS(struct('seed', seeds(si)));
    lastKey = '';
    for c = 1:nC
        sg = cfg{c,2}; alpha = 1; if sg < 0, alpha = 0; end
        key = sprintf('%g', sg);
        if ~strcmp(lastKey, key)
            g = g0; if sg > 0, g.normalSigma = sg; end
            [W, gi] = buildGraph(P_use, g);
            nComp = max(conncomp(graph(W,'upper')));
            [L, Dis] = normalizeGraph(W, alpha); [U, lam] = embedSpectral(L, 45);
            [~, ~, gg] = unique(issPartLabels(P_use(gi.idxKeep,:)));
            major = find(accumarray(gg,1) >= 60);
            lastKey = key;
        end
        o = cfg{c,4};
        switch cfg{c,5}
            case 'cheeger', Yuse = Dis * U(:, nComp+1 : nComp+15); o.phiMax = sqrt(2*lam(nComp+1));
            case 'legacy',  Yuse = Dis * U(:, 4:18);               o.maxCuts = 79;
            otherwise,      Yuse = Dis * U(:, selectModes(U, nComp, 15, cfg{c,3}));
        end
        lab = double(segmentEmbedding(Yuse, W, o));
        v  = lab > 0;
        Ct = accumarray([lab(v) gg(v)], 1, [max(lab) max(gg)]);
        PUR(c,si) = sum(max(Ct,[],2)) / sum(Ct(:));
        fr = zeros(numel(major),1);
        for q = 1:numel(major)
            col = sort(Ct(:,major(q)), 'descend'); fr(q) = find(cumsum(col) >= 0.9*sum(col), 1);
        end
        FRG(c,si) = mean(fr);
        e = evalSegmentation(uint32(lab), gg);
        ARI(c,si) = e.ARI; KK(c,si) = max(lab); COV(c,si) = mean(v);
    end
    fprintf('시드 %d 완료 (N=%d)\n', seeds(si), size(W,1));
end

fprintf('\n%-22s %-14s %6s %6s %5s %5s  %s\n', '설정', 'ARI 평균±std', '순도', '조각수', 'K', 'cov', '시드별 ARI');
for c = 1:nC
    fprintf('%-22s %.3f ± %.3f  %6.3f %6.2f %5.1f %5.2f  %s\n', cfg{c,1}, mean(ARI(c,:)), std(ARI(c,:)), ...
        mean(PUR(c,:)), mean(FRG(c,:)), mean(KK(c,:)), mean(COV(c,:)), mat2str(round(ARI(c,:),2)));
end
T = table(cfg(:,1), mean(ARI,2), std(ARI,0,2), mean(PUR,2), mean(FRG,2), mean(KK,2), mean(COV,2), ...
    'VariableNames', {'name','ARImean','ARIstd','purity','fragAvg','K','coverage'});
outDir = fullfile('out_cutting_ISS', ['seedcheck_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
names = cfg(:,1);
save(fullfile(outDir,'results.mat'), 'T','names','seeds','ARI','PUR','FRG','KK','COV');
fprintf('\n저장 완료: %s\n', outDir);
