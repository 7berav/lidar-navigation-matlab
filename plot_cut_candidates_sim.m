% run_cut_candidates_sim.m 결과를 그림으로 저장 (유용한 방식만, 방식 하나당 한 장)
%
% 입력: out_cutting_simcut/<가장 최근>_candidates/  (results.mat, labels.mat)
%       밖에서 runDir 를 주면 그 폴더, figMethods 를 주면 그 방식들만
% 출력: <runDir>/figs/<target>/<cloud>__00_GT.png         정답
%                              <cloud>__<번호>_<방식>.png 방식별 한 장
%                              <cloud>__compare.png       정답 + 방식들 (2x2)
%       <runDir>/summary_main.png  주요 방식만 고른 요약 막대그래프
%
% 기본으로 그리는 방식 (나머지는 결과가 겹치거나 탈락해서 뺐다)
%   C5_rho    : 재귀 이분할 + 병목 판정 (K 불필요)  — 현재 최선
%   K1_E1_C2  : 비율 간격 K1 + MST/NCut              — "K 를 먼저 정하는" 방식의 대표
%               (같은 K 에서 PCCA+·k-means 와 결과가 거의 같아 하나만 남김)
%   orc_E1_C2 : [참고] 정답 K + MST/NCut             — K 를 맞혔을 때의 상한
%
% 색: 예측 군집을 가장 많이 겹치는 정답 부품의 색으로 칠한다 (같은 부품에 두 번째로 걸린
%     군집부터는 다른 색). 그래서 정답 그림과 같은 색이면 맞게 자른 것이다.

if ~exist('runDir','var') || isempty(runDir)
    d = dir(fullfile('out_cutting_simcut', '*_candidates'));
    d = d(arrayfun(@(x) exist(fullfile(x.folder, x.name, 'labels.mat'), 'file') == 2, d));
    assert(~isempty(d), 'labels.mat 이 있는 *_candidates 폴더가 없습니다');
    runDir = fullfile(d(end).folder, d(end).name);
end
R = load(fullfile(runDir, 'results.mat'));
F = load(fullfile(runDir, 'labels.mat'));
T = R.T;
fprintf('입력: %s\n', runDir);

names = struct( ...
    'K1_E1_C1',  'PCCA+  (K1, 임베딩 K차원)', ...
    'K1_E1_C2',  'MST + NCut  (K1, K차원)', ...
    'K1_E3_C4',  'k-means  (K1, NJW)', ...
    'K1_E1_C6',  'MST 최장간선  (K1, K차원)', ...
    'K1_E15_C6', 'MST 최장간선  (K1, 15차원)', ...
    'C5_rho',    '재귀 이분할 + 병목 판정 (K 불필요)', ...
    'C5_K4',     '재귀 이분할 + lambda3/lambda2', ...
    'orc_E1_C2', '[참고] 정답 K + MST/NCut', ...
    'K1sc_pcca', '덩어리별 K1s + PCCA+', ...
    'K1sc_ncut', '덩어리별 K1s + MST/NCut', ...
    'K1s_pcca',  '전체 K1s + PCCA+', ...
    'orc_ncut',  '[참고] 정답 K + MST/NCut', ...
    'rec',       '재귀 이분할 + 병목 판정 (K 불필요)');
if ~exist('figMethods','var') || isempty(figMethods)
    figMethods = {'C5_rho', 'K1_E1_C2', 'orc_E1_C2'};
end
names.K1_E1_C2 = '비율 간격 K1 + MST/NCut';
pal = [lines(7); 0.55 0.55 0.55; 0.74 0.71 0.20; 0.10 0.75 0.80; 0.85 0.45 0.75; 0.45 0.25 0.10; ...
       0.20 0.20 0.55; 0.60 0.90 0.40; 1.00 0.70 0.20; 0.30 0.60 0.30];

%% 1) 점군별 그림 (onlySummary = true 를 주면 건너뜀)
if exist('onlySummary','var') && onlySummary, nFig = 0; else, nFig = numel(F.figdata); end
for i = 1:nFig
    fd = F.figdata(i);
    od = fullfile(runDir, 'figs', fd.target);
    if ~exist(od, 'dir'), mkdir(od); end
    P = double(fd.P); gt = double(fd.gt);
    [~, ~, gtl] = unique(gt);

    if isfield(fd, 'K1s'), kTxt = sprintf('K1=%d, K1s=%d', fd.K1, fd.K1s); else, kTxt = sprintf('K1=%d, K2=%d', fd.K1, fd.K2); end
    ttlGT = sprintf('%s / %s\n정답: 부품 %d개   (%s)', fd.target, fd.cloud, fd.Kgt, kTxt);
    saveOne(fullfile(od, sprintf('%s__00_GT.png', fd.cloud)), P, gtl, gtl, pal, ttlGT);

    nT = numel(figMethods) + 1; nc = ceil(sqrt(nT)); nr = ceil(nT / nc);
    fc = figure('Visible','off', 'Position',[0 0 850*nc 750*nr]);
    tl = tiledlayout(fc, nr, nc, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s / %s', fd.target, fd.cloud), 'Interpreter','none', 'FontSize', 14);
    drawOne(nexttile(tl), P, gtl, gtl, pal, sprintf('정답: 부품 %d개', fd.Kgt));

    for is = 1:numel(figMethods)
        mth = figMethods{is};
        r = T(strcmp(T.target, fd.target) & strcmp(T.cloud, fd.cloud) & strcmp(T.method, mth), :);
        lab = double(fd.labels.(mth));
        tag = '';
        if r.K > fd.Kgt, tag = sprintf('  과분할 +%d', r.K - fd.Kgt); elseif r.K < fd.Kgt, tag = sprintf('  미분할 -%d', fd.Kgt - r.K); end
        t1 = sprintf('%s\nK=%d (정답 %d)%s,  ARI=%.2f', names.(mth), r.K, fd.Kgt, tag, r.ARI);
        saveOne(fullfile(od, sprintf('%s__%02d_%s.png', fd.cloud, is, mth)), P, lab, gtl, pal, ...
                sprintf('%s / %s\n%s', fd.target, fd.cloud, t1));
        drawOne(nexttile(tl), P, lab, gtl, pal, t1);
    end
    exportgraphics(fc, fullfile(od, sprintf('%s__compare.png', fd.cloud)), 'Resolution', 160); close(fc);
    fprintf('  %s / %s\n', fd.target, fd.cloud);
end

%% 2) 주요 방식 요약 (run_cut_candidates_sim 결과일 때만)
if ~any(strcmp(T.method, 'C5_rho'))
    fprintf('저장 완료: %s\\figs\n', runDir); return
end
main = {'K1_E1_C1','K1_E1_C2','K1_E3_C4','K1_E1_C6','K1_E15_C6','E1_C3','C5_K4','C5_rho','orc_E1_C2'};
mlab = {'PCCA+','MST+NCut','k-means','최장간선(K차원)','최장간선(15차원)','HDBSCAN', ...
        '재귀 λ3/λ2','재귀 병목ρ','[참고] 정답K'};       % 한 줄로 (줄바꿈을 넣으면 눈금 이름이 밀린다)
S = zeros(numel(main), 4);
for im = 1:numel(main)
    q = strcmp(T.method, main{im}) & T.isTune == 0;
    S(im,:) = [mean(T.ARI(q)), mean(max(T.K(q)-T.Kgt(q),0)), mean(max(T.Kgt(q)-T.K(q),0)), mean(T.frag(q))];
end
f = figure('Visible','off', 'Position',[30 30 1900 1200]);
tl = tiledlayout(f, 3, 1, 'TileSpacing','compact', 'Padding','compact');
title(tl, '테스트 타깃 9종 x 점군 4개 평균  (K1 = 최대 비율 간격으로 K 선택)', 'FontSize', 13);
ax = nexttile(tl); b = bar(ax, S(:,1), 0.6); b.FaceColor = 'flat'; b.CData(end,:) = [0.7 0.7 0.7];
ylabel(ax, '평균 ARI'); ylim(ax, [0 1]); grid(ax,'on'); set(ax, 'XTick', 1:numel(main), 'XTickLabel', mlab);
text(ax, 1:numel(main), S(:,1)+0.04, compose('%.2f', S(:,1)), 'HorizontalAlignment','center');
ax = nexttile(tl); bar(ax, S(:,[2 3]), 0.6, 'stacked');
ylabel(ax, '평균 |K - 정답|'); legend(ax, {'과분할 (K 초과)','미분할 (K 부족)'}, 'Location','northwest');
grid(ax,'on'); set(ax, 'XTick', 1:numel(main), 'XTickLabel', mlab);
ax = nexttile(tl); bar(ax, S(:,4), 0.6); yline(ax, 1, '--', '정답');
ylabel(ax, '조각수 (부품 하나를 몇 조각으로)'); grid(ax,'on'); set(ax, 'XTick', 1:numel(main), 'XTickLabel', mlab);
exportgraphics(f, fullfile(runDir, 'summary_main.png'), 'Resolution', 100); close(f);
fprintf('저장 완료: %s\\figs, summary_main.png\n', runDir);

%% ===================== local functions =====================
function saveOne(file, P, lab, gtl, pal, ttl)
f = figure('Visible','off', 'Position',[0 0 1000 900]);
drawOne(axes(f), P, lab, gtl, pal, ttl);
exportgraphics(f, file, 'Resolution', 150); close(f);
end

function drawOne(ax, P, lab, gtl, pal, ttl)
% 예측 군집을 가장 많이 겹치는 정답 부품 색으로 칠함 (같은 부품의 두 번째 군집부터는 여분 색)
[~, ~, li] = unique(lab); K = max(li); G = max(gtl);
sz = accumarray(li, 1); [~, ord] = sort(sz, 'descend');
used = false(G,1); col = zeros(K,3); extra = G;
for k = ord(:)'
    g = mode(gtl(li == k));
    if ~used(g), col(k,:) = pal(mod(g-1, size(pal,1))+1, :); used(g) = true;
    else, extra = extra + 1; col(k,:) = pal(mod(extra-1, size(pal,1))+1, :); end
end
ms = 4; if numel(lab) > 12000, ms = 2.5; end
hold(ax,'on');
for k = 1:K
    m = li == k;
    scatter3(ax, P(m,1), P(m,2), P(m,3), ms, col(k,:), 'filled');
end
hold(ax,'off'); axis(ax,'equal','tight','off'); view(ax, [1 -1.2 0.8]);
title(ax, ttl, 'Interpreter','none', 'FontSize', 11);
end
