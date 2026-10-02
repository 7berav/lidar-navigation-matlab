% 합성 타깃 전체에서 "K 선택 x 임베딩 x 분할" 후보 비교
% 가중치는 C8a 와 같은 법선 가중(sigma 0.2) + alpha=1 로 고정
%
% A. K 선택
%   K1  : 최대 비율 간격  argmax lambda_{K+1}/lambda_K            (selectKRatio 'max')
%   K2  : 비율이 r0=2 를 넘는 가장 작은 K                          (selectKRatio 'first')
%   K3  : K1·K2 결과 중 사후 검증값 Upsilon = lambda_{K+1}/max_i phi(S_i) 가 큰 쪽 (사후 선택)
%   orc : 정답 부품 수 (방식 자체의 상한을 보기 위한 참고)
%   (K4 는 C5 의 정지 규칙)
% B. 임베딩
%   E1  : L_rw 앞쪽 K 개 고유벡터, 차원 정확히 K
%   E2  : diffusion map (1-lambda_i)^t psi_i,  t = 2/lambda_{K+1}.  C1 은 K 차원, C2 는 15 차원
%   E15 : 영고유벡터를 건너뛴 15 차원 (지금까지 쓰던 방식 — 차원을 줄인 효과를 보기 위한 비교)
%   E3  : L_sym 고유벡터 K 개 + 행 정규화 (NJW)
% C. 분할
%   C1 : PCCA+ (cutPCCA)                 C2 : MST 후보 + 원 그래프 NCut (cutTreeGreedy 'ncut')
%   C3 : HDBSCAN (segmentEmbedding)      C4 : k-means (n_init 20)
%   C5 : 재귀 이분할 (cutRecursive; 정지 = K4 lambda3/lambda2 > r, 또는 rho <= tau)
%   C6 : MST 최장 간선 (cutTreeGreedy 'length')
% D. 지표
%   정답 기준 : ARI, K, 순도, 조각수     V1 : K-way NCut     V2 : Upsilon, 최소 부피 비율, 최대 phi
%
% 파라미터는 튜닝 타깃(cubesat6u, cubesat6u_mid, hexsat)에서만 정함: tau=0.25, r=5. r0=2 는 사전 지정.
% 출력: out_cutting_simcut/<yyyyMMdd_HHmmss>_candidates/
%   summary.csv, spectra.mat (점군별 고유값), results.mat, labels.mat (그림용 라벨)
%   그림은 plot_cut_candidates_sim.m 으로 따로 그린다 (방식 하나당 한 장)
rng(1);

%% 0) 설정
tuneSet = {'cubesat6u', 'cubesat6u_mid', 'hexsat'};
if ~exist('targets','var') || isempty(targets)
    targets = {'cubesat6u','cubesat6u_mid','hexsat','cubesat3u_petal','spin_drum','geo_comsat', ...
               'eo_sat','telescope','soyuz','capsule','station_t','upper_stage'};
end
clouds = {'complete', 's1_mrbar_h3', 's1_oblique_h3', 's1_pvbar_h3'};
tau = 0.25; rK4 = 5; r0 = 2; kEig = 20;
shown = {'K1_E1_C1','K1_E1_C2','K1_E3_C4','K1_E1_C6','K1_E15_C6','C5_rho','C5_K4','orc_E1_C2'};   % 라벨을 저장할 방식 (그림용)

outDir = fullfile('out_cutting_simcut', [char(datetime('now','Format','yyyyMMdd_HHmmss')) '_candidates']);
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('출력: %s\n', outDir);

%% 1) 타깃 x 점군 x 방식
rows = struct([]); spectra = struct([]); figdata = struct([]);
for it = 1:numel(targets)
    d = dir(fullfile('out_cutting_sim', ['*_' targets{it}]));
    if isempty(d), fprintf('[건너뜀] %s 점군 없음\n', targets{it}); continue; end
    runDir = fullfile(d(end).folder, d(end).name);

    for ic = 1:numel(clouds)
        C = loadSimCloud(fullfile(runDir, [clouds{ic} '.mat']));
        [W, gi] = buildGraphConvex(C.P, C.sensorPos(C.view_id,:), struct('mode','unsigned','minCompSize',50));
        P = C.P(gi.idxKeep,:); gt = C.label(gi.idxKeep);
        N = size(W,1); ms = max(60, round(0.01*N));
        n0  = max(conncomp(graph(W, 'upper')));
        Kgt = numel(unique(gt));

        q  = full(sum(W,2)); Dq = spdiags(1./q, 0, N, N);
        Wt = Dq * W * Dq; Wt = (Wt + Wt.')/2;                 % alpha=1 정규화 그래프 (NCut·phi 계산용)
        [Ln, Dis] = normalizeGraph(W, 1);
        kE = min(kEig, N-2);
        [U, lam] = embedSpectral(Ln, kE);
        Psi = Dis * U;                                        % L_rw 고유벡터

        Kr = struct('K1', selectKRatio(lam, struct('rule','max')), ...
                    'K2', selectKRatio(lam, struct('rule','first','r0',r0)), 'orc', Kgt);
        spectra(end+1).target = targets{it}; spectra(end).cloud = clouds{ic}; %#ok<SAGROW>
        spectra(end).lambda = lam; spectra(end).n0 = n0; spectra(end).Kgt = Kgt; spectra(end).N = N;

        ctx = struct('target', targets{it}, 'cloud', clouds{ic}, 'isTune', double(any(strcmp(targets{it}, tuneSet))), ...
                     'N', N, 'n0', n0, 'Kgt', Kgt, 'K1', Kr.K1, 'K2', Kr.K2);
        Lb = struct();
        for kn = {'K1','K2','orc'}
            K = min(Kr.(kn{1}), kE-1);
            w = max(1 - lam, 0) .^ (2 / max(lam(min(K+1, kE)), 1e-9));           % diffusion 가중
            E1  = Psi(:, 1:K);
            E2k = Psi(:, 1:K) .* w(1:K).';
            m15 = min(15, kE - n0);
            E2m = Psi(:, 1:min(n0+m15, kE)) .* w(1:min(n0+m15, kE)).';
            E15 = Psi(:, n0+1 : n0+m15);
            E3  = U(:, 1:K); E3 = E3 ./ max(vecnorm(E3, 2, 2), eps);

            run = {
                'E1_C1',  @() cutPCCA(E1)
                'E2_C1',  @() cutPCCA(E2k)
                'E1_C2',  @() cutTreeGreedy(E1,  Wt, K, struct('score','ncut'))
                'E2_C2',  @() cutTreeGreedy(E2m, Wt, K, struct('score','ncut'))
                'E15_C2', @() cutTreeGreedy(E15, Wt, K, struct('score','ncut'))
                'E1_C6',  @() cutTreeGreedy(E1,  Wt, K, struct('score','length','minSize',ms))
                'E15_C6', @() cutTreeGreedy(E15, Wt, K, struct('score','length','minSize',ms))
                'E3_C4',  @() cutKmeans(E3, K)
            };
            for ir = 1:size(run,1)
                name = [kn{1} '_' run{ir,1}];
                [rows, Lb.(name)] = runOne(rows, ctx, name, kn{1}, run{ir,1}, run{ir,2}, gt, Wt, lam);
            end
        end
        [rows, Lb.C5_K4]  = runOne(rows, ctx, 'C5_K4',  'K4',  'C5', @() cutRecursive(W, struct('stop','ratio','r',rK4,'minSize',ms)), gt, Wt, lam);
        [rows, Lb.C5_rho] = runOne(rows, ctx, 'C5_rho', 'rho', 'C5', @() cutRecursive(W, struct('tau',tau,'minSize',ms)), gt, Wt, lam);
        Eh = Psi(:, 1:min(Kr.K1, kE-1));
        [rows, Lb.E1_C3]  = runOne(rows, ctx, 'E1_C3',  'hdb', 'E1_C3', @() segmentEmbedding(Eh, W, ...
                                   struct('ky',min(40,N-2),'minSize',ms,'cutMethod','hdbscan','mreachK',15)), gt, Wt, lam);

        fd = struct('target', targets{it}, 'cloud', clouds{ic}, 'P', single(P), 'gt', uint8(gt), ...
                    'Kgt', Kgt, 'K1', Kr.K1, 'K2', Kr.K2, 'labels', struct());
        here = rows(strcmp({rows.target}, targets{it}) & strcmp({rows.cloud}, clouds{ic}));
        fprintf('%-16s %-14s Kgt=%2d K1=%2d K2=%2d |', targets{it}, clouds{ic}, Kgt, Kr.K1, Kr.K2);
        for is = 1:numel(shown)
            r = here(strcmp({here.method}, shown{is}));
            fd.labels.(shown{is}) = uint8(Lb.(shown{is}));
            fprintf(' %s %d/%.2f', shown{is}, r.K, r.ARI);
        end
        figdata = [figdata; fd]; %#ok<AGROW>
        fprintf('\n');
    end
    writetable(struct2table(rows), fullfile(outDir, 'summary.csv'));
    save(fullfile(outDir, 'spectra.mat'), 'spectra');
    save(fullfile(outDir, 'labels.mat'), 'figdata', 'shown');
end

%% 2) K3 (사후 검증): 같은 방식의 K1·K2 결과 중 Upsilon 이 큰 쪽
T = struct2table(rows);
k1 = T(strcmp(T.kRule,'K1'), :); add = k1([],:);
for i = 1:height(k1)
    j = find(strcmp(T.kRule,'K2') & strcmp(T.target,k1.target{i}) & strcmp(T.cloud,k1.cloud{i}) & strcmp(T.base,k1.base{i}));
    pick = k1(i,:);
    if ~isempty(j) && T.upsilon(j) > pick.upsilon, pick = T(j,:); end
    pick.kRule = {'K3'}; pick.method = {['K3_' pick.base{1}]};
    add = [add; pick]; %#ok<AGROW>
end
T = [T; add];
writetable(T, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'results.mat'), 'T', 'spectra', 'tau', 'rK4', 'r0', 'tuneSet', 'targets', 'clouds');

%% 3) 요약
fprintf('\n[K 선택 정확도] 48개 점군 (정답 = 보이는 부품 수)\n');
Kg = [spectra.Kgt].';
kr = struct('name', {'K1 최대비율','K2 r0=1.5','K2 r0=2','K2 r0=3','차이간격(기존)'}, 'K', cell(1,5));
for i = 1:numel(spectra)
    l = spectra(i).lambda; nz = spectra(i).n0;
    kr(1).K(i,1) = selectKRatio(l, struct('rule','max'));
    kr(2).K(i,1) = selectKRatio(l, struct('rule','first','r0',1.5));
    kr(3).K(i,1) = selectKRatio(l, struct('rule','first','r0',2));
    kr(4).K(i,1) = selectKRatio(l, struct('rule','first','r0',3));
    g = diff(l); kk = max(nz,1):min(15, numel(g)); [~, ig] = max(g(kk)); kr(5).K(i,1) = kk(ig);
end
for i = 1:numel(kr)
    fprintf('  %-16s 정확 %2d/%d, 과대 평균 %.2f (최대 %d), 과소 평균 %.2f (최대 %d)\n', kr(i).name, nnz(kr(i).K == Kg), numel(Kg), ...
        mean(max(kr(i).K-Kg,0)), max(kr(i).K-Kg), mean(max(Kg-kr(i).K,0)), max(Kg-kr(i).K));
end

methods = unique(T.method, 'stable');
Sm = zeros(numel(methods), 10);
fprintf('\n%-14s | 튜닝: ARI  과분할 미분할 | 테스트: ARI  과분할 미분할 조각수  순도   NCut  Upsilon\n', '방식');
for im = 1:numel(methods)
    a = strcmp(T.method, methods{im}) & T.isTune == 1;
    b = strcmp(T.method, methods{im}) & T.isTune == 0;
    Sm(im,:) = [mean(T.ARI(a)), mean(max(T.K(a)-T.Kgt(a),0)), mean(max(T.Kgt(a)-T.K(a),0)), ...
                mean(T.ARI(b)), mean(max(T.K(b)-T.Kgt(b),0)), mean(max(T.Kgt(b)-T.K(b),0)), ...
                mean(T.frag(b)), mean(T.purity(b)), mean(T.ncutK(b)), median(T.upsilon(b), 'omitnan')];
    fprintf('%-14s |      %.3f  %5.2f  %5.2f |        %.3f  %5.2f  %5.2f  %5.2f  %.3f  %.4f  %6.2f\n', methods{im}, Sm(im,:));
end

fprintf('\n저장 완료: %s  (그림은 plot_cut_candidates_sim.m)\n', outDir);

%% ===================== local functions =====================
function [rows, lab] = runOne(rows, ctx, name, kRule, base, fn, gt, Wt, lam)
t = tic;
try
    lab = fn();
catch err
    warning('%s 실패: %s', name, err.message);
    lab = zeros(numel(gt), 1, 'uint32');
end
tm = toc(t);
lab = double(lab);
if any(lab == 0), lab(lab == 0) = max(lab) + 1; end        % 미할당은 별도 군집으로
[~, ~, li] = unique(lab); [gu, ~, gidx] = unique(double(gt));
K  = max(li);
e  = evalSegmentation(uint32(li), gt);
Ct = accumarray([li gidx], 1);
fr = zeros(numel(gu), 1);
for qn = 1:numel(gu)
    col = sort(Ct(:,qn), 'descend'); fr(qn) = find(cumsum(col) >= 0.9*sum(col), 1);
end
major = sum(Ct, 1).' >= 30;

d = full(sum(Wt,2)); volTot = sum(d);
vol = accumarray(li, d, [K 1]);
[ei, ej, ew] = find(Wt);
assoc = accumarray(li(ei(li(ei) == li(ej))), ew(li(ei) == li(ej)), [K 1]);
cutv  = vol - assoc;
phi   = cutv ./ max(min(vol, volTot - vol), eps);
ups   = NaN;
if K >= 2 && K + 1 <= numel(lam), ups = lam(K+1) / max(max(phi), eps); end

r = ctx;
r.method = name; r.kRule = kRule; r.base = base;
r.K = K; r.ARI = e.ARI; r.purity = sum(max(Ct, [], 2)) / sum(Ct(:)); r.frag = mean(fr(major));
r.ncutK = sum(cutv ./ max(vol, eps)); r.upsilon = ups; r.minVolFrac = min(vol) / volTot; r.maxPhi = max(phi);
r.timeSec = tm;
rows = [rows; r];
lab = uint32(li);
end

function lab = cutKmeans(X, K)
if K <= 1, lab = ones(size(X,1), 1, 'uint32'); return; end
lab = uint32(kmeans(X, K, 'Replicates', 20, 'MaxIter', 300));
end
