% 점군별 진단 그림: 스펙트럼과 K 옵션들 + 앞쪽 고유벡터 3개(psi2..psi4) 공간의 MST
%
% 한 장에 4칸
%   (1) 고유값 lambda_k (로그 축) 와 비율 lambda_{k+1}/lambda_k (막대)
%       세로선: 정답 K, K1(최대 비율), K1s(거의 0 건너뛰기), 차이 간격(기존)
%       가로선: "거의 0" 기준 (nearRel * lambda_16), 비율 문턱 rMin = 2
%   (2) psi2, psi3, psi4 좌표에 점(색 = 정답 부품) + MST(회색) + 잘린 간선(빨강, 굵게)
%       MST 는 이 3개 벡터로만 만들고(임베딩 kNN + 원 그래프 간선의 최소신장 트리),
%       원 그래프 NCut 이 가장 적게 늘어나는 간선부터 K1s-1 번 자른다 (cutTreeGreedy 'ncut')
%   (3) 같은 MST 와 잘린 간선을 3D 실제 좌표에 그린 것
%   (4) 그 결과 (색 = 가장 많이 겹치는 정답 부품)
%
% 가중치: C8a 와 같은 법선 가중(sigma 0.2) + alpha=1 (run_cut_oneshot_sim 과 동일)
% 출력: out_cutting_simcut/<yyyyMMdd_HHmmss>_spectrum/<target>__<cloud>.png
%       밖에서 targets / clouds 를 주면 그것만

if ~exist('targets','var') || isempty(targets)
    targets = {'cubesat6u','cubesat6u_mid','hexsat','cubesat3u_petal','spin_drum','geo_comsat', ...
               'eo_sat','telescope','soyuz','capsule','station_t','upper_stage'};
end
if ~exist('clouds','var') || isempty(clouds)
    clouds = {'complete', 's1_mrbar_h3', 's1_oblique_h3', 's1_pvbar_h3'};
end
kEig = 20; nearRel = 0.01; rMin = 2;
outDir = fullfile('out_cutting_simcut', [char(datetime('now','Format','yyyyMMdd_HHmmss')) '_spectrum']);
if ~exist(outDir, 'dir'), mkdir(outDir); end
pal = [lines(7); 0.55 0.55 0.55; 0.74 0.71 0.20; 0.10 0.75 0.80; 0.85 0.45 0.75; 0.45 0.25 0.10];

for it = 1:numel(targets)
    d = dir(fullfile('out_cutting_sim', ['*_' targets{it}]));
    if isempty(d), continue; end
    runDir = fullfile(d(end).folder, d(end).name);
    for ic = 1:numel(clouds)
        C = loadSimCloud(fullfile(runDir, [clouds{ic} '.mat']));
        [W, gi] = buildGraphConvex(C.P, C.sensorPos(C.view_id,:), struct('mode','unsigned','minCompSize',50));
        P = C.P(gi.idxKeep,:); gt = double(C.label(gi.idxKeep)); [~, ~, gtl] = unique(gt);
        N = size(W,1); n0 = max(conncomp(graph(W,'upper'))); Kgt = max(gtl);
        q = full(sum(W,2)); Dq = spdiags(1./q, 0, N, N); Wt = Dq*W*Dq; Wt = (Wt+Wt.')/2;
        [Ln, Dis] = normalizeGraph(W, 1);
        kE = min(kEig, N-2);
        [U, lam] = embedSpectral(Ln, kE);
        Psi = Dis * U;

        K1  = selectKRatio(lam, struct('rule','max'));
        [K1s, inf1s] = selectKRatio(lam, struct('rule','skipnear','nearRel',nearRel,'rMin',rMin));
        g = diff(lam); kk = max(n0,1):min(15, numel(g)); [~, ig] = max(g(kk)); Kabs = kk(ig);
        ratio = lam(2:end) ./ max(lam(1:end-1), 1e-12);
        nearThr = nearRel * lam(min(16, kE));

        Y3 = Psi(:, 2:4);
        [lab, info] = cutTreeGreedy(Y3, Wt, K1s, struct('score','ncut'));
        e = evalSegmentation(lab, uint32(gt));

        f = figure('Visible','off', 'Position',[0 0 1900 1500]);
        tl = tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');
        title(tl, sprintf('%s / %s   (N=%d, 끊어진 덩어리 %d, 정답 부품 %d)', targets{it}, clouds{ic}, N, n0, Kgt), ...
              'Interpreter','none', 'FontSize', 14);

        % (1) 스펙트럼
        ax = nexttile(tl);
        % (로그 축 + 반투명 막대 조합은 보이지 않는 그림에서 렌더링되지 않아 log10 값을 선형 축에 그림)
        yyaxis(ax, 'left');
        plot(ax, 1:kE, log10(max(lam, 1e-12)), 'o-', 'LineWidth', 1.4, 'MarkerFaceColor', [0 0.45 0.74]); hold(ax, 'on');
        yline(ax, log10(nearThr), ':', '거의 0 기준', 'LabelHorizontalAlignment','left', 'LineWidth', 1.2);
        ylabel(ax, 'log_{10} \lambda_k');
        yyaxis(ax, 'right');
        bar(ax, (1:kE-1) + 0.5, ratio, 0.5, 'FaceColor', [0.98 0.80 0.70], 'EdgeColor', 'none');
        yline(ax, rMin, '--', sprintf('비율 문턱 %g', rMin), 'LabelHorizontalAlignment','left');
        ylabel(ax, '\lambda_{k+1} / \lambda_k'); ylim(ax, [0 min(max(ratio(2:end))*1.15+1, 30)]);
        hold(ax, 'off');
        mk = {Kgt, '정답', [0 0.6 0], '-'; K1, 'K1', [0.85 0.33 0.1], '--'; K1s, 'K1s', [0.49 0.18 0.56], '-.'; Kabs, '차이간격', [0.4 0.4 0.4], ':'};
        for m = 1:size(mk,1)
            xline(ax, mk{m,1} + 0.5, mk{m,4}, sprintf('%s=%d', mk{m,2}, mk{m,1}), 'Color', mk{m,3}, 'LineWidth', 1.6, ...
                  'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        end
        grid(ax, 'on'); xlim(ax, [0.5 kE+0.5]); xlabel(ax, 'k');
        title(ax, sprintf('스펙트럼: 거의 0 인 고유값 %d개 (진짜 0 = 덩어리 %d개)', inf1s.nNear, n0));

        % (2) psi2..psi4 공간의 MST
        Ea = info.E(info.alive, :); Ec = info.E(~info.alive, :);
        ax = nexttile(tl); hold(ax, 'on');
        drawEdges(ax, Y3, Ea, [0.6 0.6 0.6], 0.5);
        cols = pal(mod(gtl-1, size(pal,1))+1, :);
        scatter3(ax, Y3(:,1), Y3(:,2), Y3(:,3), 6, cols, 'filled');
        drawEdges(ax, Y3, Ec, [1 0 0], 3);
        if ~isempty(Ec)
            mid = (Y3(Ec(:,1),:) + Y3(Ec(:,2),:)) / 2;
            scatter3(ax, mid(:,1), mid(:,2), mid(:,3), 90, 'r', 'LineWidth', 1.5);
        end
        hold(ax, 'off'); grid(ax, 'on'); view(ax, 3); axis(ax, 'tight');
        xlabel(ax, '\psi_2'); ylabel(ax, '\psi_3'); zlabel(ax, '\psi_4');
        title(ax, sprintf('고유벡터 \\psi_2~\\psi_4 공간의 MST (회색), 잘린 간선 %d개 (빨강 원)', size(Ec,1)));

        % (3) 같은 MST 를 실제 좌표에
        ax = nexttile(tl); hold(ax, 'on');
        drawEdges(ax, P, Ea, [0.55 0.55 0.55], 0.4);
        drawEdges(ax, P, Ec, [1 0 0], 3);
        if ~isempty(Ec)
            mid = (P(Ec(:,1),:) + P(Ec(:,2),:)) / 2;
            scatter3(ax, mid(:,1), mid(:,2), mid(:,3), 120, 'r', 'LineWidth', 1.5);
        end
        hold(ax, 'off'); axis(ax, 'equal', 'tight'); grid(ax, 'on'); view(ax, [1 -1.2 0.8]);
        title(ax, '같은 MST 를 실제 3D 좌표에 (빨강 원 = 잘린 곳)');

        % (4) 결과
        ax = nexttile(tl);
        drawLab(ax, P, double(lab), gtl, pal);
        title(ax, sprintf('결과: MST+NCut, K = K1s = %d (정답 %d), ARI = %.2f', K1s, Kgt, e.ARI));

        fn = fullfile(outDir, sprintf('%s__%s.png', targets{it}, clouds{ic}));
        exportgraphics(f, fn, 'Resolution', 110); close(f);
        fprintf('%-16s %-14s 정답 %2d | K1 %2d K1s %2d 차이간격 %2d | 거의0 %d (덩어리 %d) | ARI %.2f\n', ...
            targets{it}, clouds{ic}, Kgt, K1, K1s, Kabs, inf1s.nNear, n0, e.ARI);
    end
end
fprintf('저장 완료: %s\n', outDir);

function drawEdges(ax, X, E, col, lw)
if isempty(E), return; end
n = size(E,1);
xs = [X(E(:,1),1), X(E(:,2),1), nan(n,1)].';
ys = [X(E(:,1),2), X(E(:,2),2), nan(n,1)].';
zs = [X(E(:,1),3), X(E(:,2),3), nan(n,1)].';
plot3(ax, xs(:), ys(:), zs(:), '-', 'Color', col, 'LineWidth', lw);
end

function drawLab(ax, P, lab, gtl, pal)
[~, ~, li] = unique(lab); K = max(li); G = max(gtl);
sz = accumarray(li, 1); [~, ord] = sort(sz, 'descend');
used = false(G,1); col = zeros(K,3); extra = G;
for k = ord(:)'
    g = mode(gtl(li == k));
    if ~used(g), col(k,:) = pal(mod(g-1, size(pal,1))+1, :); used(g) = true;
    else, extra = extra + 1; col(k,:) = pal(mod(extra-1, size(pal,1))+1, :); end
end
hold(ax, 'on');
for k = 1:K
    m = li == k;
    scatter3(ax, P(m,1), P(m,2), P(m,3), 4, col(k,:), 'filled');
end
hold(ax, 'off'); axis(ax, 'equal', 'tight', 'off'); view(ax, [1 -1.2 0.8]);
end
