syms x y z
order  = 2;
TermsC = homogeneFischerTerms(order);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
nT = length(TermsC);
%N2 = length(TermsC);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);

rot_each   = true; 
range  = [-60 60 -35 40 -30 28];
range_local  = [-10 10 -7 7 -10 10];
faceColor = [0.95, 0.72, 0.5];
faceAlpha = 0.50; meshDen_Data = 80;
warning('off','MATLAB:nearlySingularMatrix')

%%
% === 5) 결과 버퍼 ===
%Res = [];   % 테이블로 모을 예정
%row = 0;
%%
%sigma_list = [ 0.05 ];     % 가우시안 노이즈 표준편차
%q_list     = [0.70 0.5 0.0 0.1 0.30];  
%%
%{
% === 6) 실험 루프 ===
for s = sigma_list
  % 노이즈 추가
  Pnoisy = PP_use + randn(size(PP_use))*s;

  for q = q_list
    P = Pnoisy;



    % 가시율: z-축 프로젝션 상위 q 분위수만 유지
    zproj = P(:,3);
    tau   = quantile(zproj, q);    % 또는 prctile(zproj, q*100)
    keep  = zproj >= tau;
    P_use = P(keep, :);

    if size(P_use,1) < nT*2
        fprintf('[skip] sigma=%.3f, q=%.2f → 표본 너무 적음(N=%d)\n', s, q, size(P_use,1));
        continue;
    end

    % 단일 피팅 실행 (회귀1회 + 로컬1회)
    [Disp, Beta, inlierMask, Log] = ...
        PoliNavigationSolver3_Fit(P_use, order, nT, Funcs1, Grads1, fitPar);
    
    % Log = [preScore, postScore, raw_ms, local_ms]
    % 결과 집계
    row = row + 1;
    Res(row).sigma     = s;

    Res(row).q         = q;
    Res(row).iter      = fitPar.locIters;    
    Res(row).N         = size(P_use,1);
    Res(row).preScore  = Log(1);
    Res(row).postScore = Log(2);
    Res(row).t_raw_ms  = Log(3);
    Res(row).t_loc_ms  = Log(4);
    Res(row).inlierR   = mean(inlierMask);   % 참고치
    
    fprintf('σ=%.3f, q=%.2f, N=%5d | pre=%.4f → post=%.4f | raw=%.1fms, loc=%.1fms\n', ...
        s, q, Res(row).N, Res(row).preScore, Res(row).postScore, Res(row).t_raw_ms, Res(row).t_loc_ms);
  
    % 입력: P_use (Nx3), Beta (nT×1), ResultDisp(1×3), Funcs, Grads
    P_result = P - Disp;                          % 모델 좌표계로 평행이동
    phi   = Funcs1([P_result(:,1), P_result(:,2), P_result(:,3)]);   % N×nT
    Fval  = phi * Beta - 1;                              % N×1, 암시적 잔차
    
    Gstack = Grads1(P_result);                               % (3N)×nT
    GradVal = Gstack * Beta;                             % (3N)×1
    Nn = size(P_result,1);
    gn = sqrt( max(1e-12, ...
          GradVal(1:Nn).^2 + GradVal(Nn+1:2*Nn).^2 + GradVal(2*Nn+1:3*Nn).^2) );
    
    d_sdf = Fval ./ gn;                % 부호 거리
    %d_sdf = abs(Fval) ./ gn;                % 절대 거리

    SDF_RMS = sqrt(mean(d_sdf.^2));
    SDF_P95 = prctile(d_sdf, 95);
    Res(row).Dispx     = Disp(1);
    Res(row).Dispy     = Disp(2);
    Res(row).Dispz     = Disp(3);
    Res(row).SDF_RMS   = SDF_RMS;
    Res(row).SDF_P95   = SDF_P95;
  end
end

% === 7) 표로 정리 ===
T = struct2table(Res);
%}
%%
%% 0) 데이터셋/회전/리스트 정의
% --- 점군 여러 개 준비(원하는 생성기로 교체 가능) ---
P_sets = {
  generateRandomPointsOnHexagonPrism(10000);
  generateRandomPointsOnHexagonPrism(10000);

  generateRandomPointsOnCylinder(10000);
  generateRandomPointsOnCylinder(10000);
  generateRandomPointsOnCube(10000);
  generateRandomPointsOnCube(10000);
  generateRandomPointsOnSphere(10000);
  generateRandomPointsOnEllipsoid(10000,0.8,1.0,1.2);
  generateRandomPointsOnCube(10000);
  generateRandomPointsOnHexagonPrism(10000);
};
nP = numel(P_sets);

% 점 수정
for pid = 2:2:6
    P = P_sets{pid};
    scale = [0.8, 1.0, 1.2];          % x 0.8배, y 1배, z 2배
    P = P .* scale;   
    P_sets{pid} = P;
end
for pid = 9:10
    P = P_sets{pid};
    scale = [0.6, 0.6, 1.5];          % x 0.8배, y 1배, z 2배
    P = P .* scale;   
    P_sets{pid} = P;
end


% --- 쿼터니언(회전) 리스트: 10개 샘플 (Fibonacci-유사 방향 -> 회전축 + 고정 각) ---
nQ = 10;
q  = randn(nQ, 4);
q  = q ./ vecnorm(q, 2, 2);      % 각 행 단위화

Q_list = cell(nQ,1);
for qi = 1:nQ
    w = q(qi,1); x = q(qi,2); y = q(qi,3); z = q(qi,4);
    R = [1-2*(y^2+z^2)  2*(x*y - z*w)  2*(x*z + y*w);
         2*(x*y + z*w)  1-2*(x^2+z^2)  2*(y*z - x*w);
         2*(x*z - y*w)  2*(y*z + x*w)  1-2*(x^2+y^2)];
    Q_list{qi} = R;
end
% --- 실험 리스트 ---
sigma_list = [0.02 0.05];                % 노이즈 σ
q_list     = [0.40 0.33 0.27 0.10 0.15 0.2 0.0];  % 상위 q 분위수 컷(남는 비율 ≈ 1-q)
it_list    = [4];
% --- 피팅/스코어 파라미터 ---
fitPar = struct( ...
  'scoreType' ,'similarity', ...
  'conf', 0.96, ...
  'thresh', 0.40, ...       % 요청대로 0.4 고정
  'locIters', 4, ...
  'minInlierRatio', 0.65, ...
  'momentum', 0.66, ...
  'damping', 0.85, ...
  'reg', 1e-6, ...
  'sigma', 0.13);           % similarity 스코어 폭

%% 2) 결과 버퍼
Res = []; row = 0;

%% 3) 전체 중첩 루프: 데이터셋 × 회전 × σ × q
for it = it_list 
    fitPar.locIters = it;
    for pid = 1:nP
      P0 = P_sets{pid};
    
      for Qid = 1:nQ
        R = Q_list{Qid};
        Prot = (R * P0.').';                           % 회전 적용
    
        for s = sigma_list
          Pnoisy = Prot + randn(size(Prot))*s;         % σ-노이즈
    
          for q = q_list
              
                % 가시율: z-축 상위 q 분위수만 유지
                zproj = Pnoisy(:,3);
                tau   = quantile(zproj, q);
                keep  = (zproj >= tau);
                P_use = Pnoisy(keep, :);
        
                if size(P_use,1) < nT*2
                    fprintf('[skip] pid=%d, Qid=%d, sigma=%.3f, q=%.2f → N=%d\n', ...
                            pid, Qid, s, q, size(P_use,1));
                    continue;
                end
        
                % 단일 피팅 (회귀1회 + 로컬1회)
                [Disp, Beta, inlierMask, Log] = ...
                    PoliNavigationSolver3_Fit(P_use, order, nT, Funcs1, Grads1, fitPar);
        
                % 외부 SDF (signed)
                P_result = Pnoisy - Disp;   % 모델 좌표계
                Phi      = Funcs1([P_result(:,1), P_result(:,2), P_result(:,3)]);  % N×nT
                Fval     = Phi * Beta - 1;                                       % N×1
                Gstack   = Grads1(P_result);                                     % (3N)×nT
                gv       = Gstack * Beta;  Nn = size(P_result,1);
                gn       = sqrt( max(1e-12, gv(1:Nn).^2 + gv(Nn+1:2*Nn).^2 + gv(2*Nn+1:3*Nn).^2) );
                d_sdf    = Fval ./ gn;
                
                inlierMask_all = abs(Fval) <= fitPar.thresh;
                SDF_RMS  = sqrt(mean(d_sdf.^2));
                SDF_P50  = prctile(d_sdf, 50);
        
                % 결과 저장 (필수 요약만)
                row = row + 1;
                %Res(row).it        = it;
                Res(row).pid       = pid;
                Res(row).Qid       = Qid;
                Res(row).sigma     = s;
                Res(row).q         = q;
                Res(row).iter      = fitPar.locIters;
                Res(row).N         = size(P_use,1);
                Res(row).preScore  = Log(1);
                Res(row).postScore = Log(2);
                Res(row).t_raw_ms  = Log(3);
                Res(row).t_loc_ms  = Log(4);
                Res(row).inlierR   = mean(inlierMask_all);
                Res(row).Dispx     = Disp(1);
                Res(row).Dispy     = Disp(2);
                Res(row).Dispz     = Disp(3);
                Res(row).SDF_RMS   = SDF_RMS;
                Res(row).SDF_P50   = SDF_P50;
        
                fprintf('pid=%d Q=%2d σ=%.3f q=%.2f N=%5d | post=%.4f | raw=%.1fms loc=%.1fms\n', ...
                  pid, Qid, s, q, Res(row).N, Log(2), Log(3), Log(4));
              
          end
        end
      end
    end
end

T = struct2table(Res);
%%
d = d_sdf;

figure(21); clf;
histogram(d, 'NumBins', 80);    % 필요하면 'Normalization','pdf'
grid on; box on;
xlabel('d\_sdf'); ylabel('count');
xlim([-1 1]);
%title('Histogram of d\_sdf');
%%
P_use_unbias =  P_use - Disp;
P_total_unbias = Pnoisy - Disp;
if ~isempty(Beta)
    f1_RAN = TermsC * Beta;
    Funcs1_RAN = matlabFunction(f1_RAN);
    
    P_inlier = P_use_unbias(logical(inlierMask),:);
    %PP_inlier = P_use(logical(inlierMask), :);
    P_inlier_unbias = P_inlier - Disp;
    Value_inlier_unbias = Funcs1_RAN(P_inlier_unbias(:,1),P_inlier_unbias(:,2),P_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(P_inlier_unbias(:,1),P_inlier_unbias(:,2),P_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,[-20 20 -10 10 -30 30],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 150);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, -0.2]);
    axis equal
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

    
    figure(3)
    scatter3(P_total_unbias(:,1),P_total_unbias(:,2),P_total_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(P_inlier_unbias(:,1),P_inlier_unbias(:,2),P_inlier_unbias(:,3),3,'b','filled');
    
    fimplicit3(f1_RAN-1,[-20 20 -20 20 -30 30],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor','none','FaceAlpha',0.5, 'MeshDensity', 150);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, -0.2]);
    axis equal
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

end

%% 1-2. q–σ에 따른 SDF
it_ref = 4;        % 분석에 사용할 iter (가장 큰 값)
T1 = T(T.iter == it_ref, :); % 해당 iter만 사용

G = groupsummary(T1, {'sigma','q'}, {'mean','median'}, 'SDF_RMS'); 
% G.mean_SDF_RMS : 같은 (sigma, q) 조합에 대한 SDF_RMS 평균
q_max_main = 0.6;

sigmas = unique(G.sigma);
figure; hold on; grid on;
for i = 1:numel(sigmas)
    s  = sigmas(i);
    idx = (G.sigma == s)& (G.q <= q_max_main);
    qv  = G.q(idx);
    yv  = G.median_SDF_RMS(idx);
    
    % q 오름차순 정렬
    [qv, order] = sort(qv);
    yv = yv(order);
    
    plot(qv, yv, '-o', 'DisplayName', sprintf('\\sigma = %.3f', s));
end
xlabel('visibility ');
ylabel('SDF ');
%title(sprintf('SDF\\_RMS vs q (iter = %d)', it_ref));
legend('Location','best');
hold off;


%% ==== 1-2. q–σ에 따른 inlierR 선 그래프 ====
G2 = groupsummary(T1, {'sigma','q'}, 'mean', 'inlierR'); 
% G2.mean_inlierR : 같은 (sigma, q) 조합의 inlier 비율 평균

sigmas = unique(G2.sigma);
figure; hold on; grid on;
for i = 1:numel(sigmas)
    s  = sigmas(i);
    idx = (G2.sigma == s)& (G2.q <= q_max_main);
    qv  = G2.q(idx);
    yv  = G2.mean_inlierR(idx);
    
    [qv, order] = sort(qv);
    yv = yv(order);
    
    plot(qv, yv, '-o', 'DisplayName', sprintf('\\sigma = %.3f', s));
end
xlabel('visibility');
ylabel('inlier ratio');
%title(sprintf('inlierR vs q (iter = %d)', it_ref));
ylim([0.4 1.2])
legend('Location','best');
hold off;

%% ==== 2) iter–q에 따른 SDF_RMS, Dispz ====

sigma_ref = min(T.sigma);   % 보고 싶은 sigma 선택 (예: 가장 작은 sigma)
T2 = T(abs(T.sigma - sigma_ref) < 1e-12, :);  % 해당 sigma만 필터링
q_max_main = 0.6;
% 2-1. iter–q에 따른 SDF_RMS
G_iter_sdf = groupsummary(T2, {'iter','q'}, {'mean','median'}, 'SDF_RMS');
% G_iter_sdf.mean_SDF_RMS : 같은 (iter, q) 조합에서의 SDF_RMS 평균

q_list = unique(G_iter_sdf.q);
figure; hold on; grid on;
for i = 1:numel(q_list)
    qv = q_list(i);
    idx = (G_iter_sdf.q == qv)& (G_iter_sdf.q <= q_max_main);
    itv = G_iter_sdf.iter(idx);
    yv  = G_iter_sdf.median_SDF_RMS(idx);
    
    [itv, order] = sort(itv);
    yv = yv(order);
    
    plot(itv, yv, '-o', 'DisplayName', sprintf('q = %.2f', qv));
end
xlabel('number of local iterations (iter)');
ylabel('SDF');
%title(sprintf('SDF\\_RMS vs iter (\\sigma = %.3f)', sigma_ref));
legend('Location','best');
xlim([ 1 8]); 
hold off;

% 2-2. iter–q에 따른 Dispz
G_iter_dispz = groupsummary(T2, {'iter','q'}, {'mean','median'}, 'Dispz');
% G_iter_dispz.mean_Dispz : 같은 (iter, q) 조합에서의 z-방향 변위 평균

q_list = unique(G_iter_dispz.q);
figure; hold on; grid on;
for i = 1:numel(q_list)
    qv = q_list(i);
    idx = (G_iter_dispz.q == qv)& (G_iter_dispz.q <= q_max_main);
    itv = G_iter_dispz.iter(idx);
    yv  = G_iter_dispz.median_Dispz(idx);
    
    [itv, order] = sort(itv);
    yv = yv(order);
    
    plot(itv, yv, '-o', 'DisplayName', sprintf('q = %.2f', qv));
end
xlabel('number of local iterations ');
ylabel('z coordinate of center');
%title(sprintf('Dispz vs iter (\\sigma = %.3f)', sigma_ref));
legend('Location','best');
xlim([ 1 8]); %ylim([-0.4 1.2]);
hold off;
%saveas(gcf, 'monte_fitting\ord2\graph_zcenter_iter_3', 'svg');
%savefig(gcf, 'monte_fitting\ord2\hist_sdf.fig');
%exportgraphics(gcf, 'monte_fitting\ord2\graph_zcenter_iter_3_2.pdf', 'ContentType', 'vector');


% 2-23. iter–q에 따른 inlier
G_iter_inlier = groupsummary(T2, {'iter','q'}, {'mean','median'}, 'inlierR');

q_list = unique(G_iter_inlier.q);
figure; hold on; grid on;
for i = 1:numel(q_list)
    qv = q_list(i);
    idx = (G_iter_inlier.q == qv)& (G_iter_inlier.q <= q_max_main);
    itv = G_iter_inlier.iter(idx);
    yv  = G_iter_inlier.median_inlierR(idx);
    
    [itv, order] = sort(itv);
    yv = yv(order);
    
    plot(itv, yv, '-o', 'DisplayName', sprintf('q = %.2f', qv));
end
xlabel('number of local iterations ');
ylabel('inlier ratio');
%title(sprintf('Dispz vs iter (\\sigma = %.3f)', sigma_ref));
legend('Location','best');
xlim([ 1 8]); %ylim([-0.4 1.2]);
hold off;

%% 1-3. 시간 통계 (히스토그램 + 가시율별 통계)

% 시간 분석에 사용할 iter 선택 (위 SDF 분석과 동일하게 사용)
it_ref_time = 4;  % 필요시 2, 7 등으로 변경 가능
T_time = T(T.iter == it_ref_time, :);   % 해당 iter만 사용
% 전체 iter에 대해 보고 싶으면:
% T_time = T;

% ----- Fig6: 전체 시간 히스토그램(겹쳐, 선형축) -----
raw_ms = T_time.t_raw_ms;
loc_ms = T_time.t_loc_ms;

% NaN, Inf, 음수 제거
raw_ms = raw_ms(isfinite(raw_ms) & raw_ms > 0);
loc_ms = loc_ms(isfinite(loc_ms) & loc_ms > 0);

has_time_raw = ~isempty(raw_ms);
has_time_loc = ~isempty(loc_ms);

if has_time_raw || has_time_loc
    figure(6); clf; hold on; grid on;

    if has_time_raw
        % FD(Freedman–Diaconis) 규칙 기반 bin 경계 계산
        [~, edges_raw] = histcounts(raw_ms, 'BinMethod', 'fd');
        histogram(raw_ms, edges_raw, ...
            'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, ...
            'DisplayName', 'raw (1st regression)');
    end

    if has_time_loc
        [~, edges_loc] = histcounts(loc_ms, 'BinMethod', 'fd');
        histogram(loc_ms, edges_loc, ...
            'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, ...
            'DisplayName', 'local (local iters)');
    end

    xlabel('time [ms]');
    ylabel('count');
    %title(sprintf('Iteration time (iter = %d)', it_ref_time));

    % 상위 99 percentile 정도까지로 x축 제한 (긴 꼬리 잘라서 보기 좋게)
    all_t = [raw_ms; loc_ms];
    if ~isempty(all_t)
        t_max = prctile(all_t, 99);   % 상위 1%는 잘라냄
        xlim([0, t_max]);
    end

    legend('Location','best');
    hold off;
end



%% q별 시간 그래프
q_max_time = 0.6;   % 상위 q만 보고 싶으면 제한 (예: 0.6까지)

% sigma, q별로 t_raw_ms, t_loc_ms에 대한 mean/median 계산
G_time = groupsummary(T_time, {'sigma','q'}, ...
    {'mean','median'}, {'t_raw_ms','t_loc_ms'});
% G_time에는 다음 변수들이 생김:
%   G_time.mean_t_raw_ms,  G_time.median_t_raw_ms
%   G_time.mean_t_loc_ms,  G_time.median_t_loc_ms

sigmas = unique(G_time.sigma);

figure; hold on; grid on;
for i = 1:numel(sigmas)
    s  = sigmas(i);
    idx = (G_time.sigma == s) & (G_time.q <= q_max_time);

    qv = G_time.q(idx);
    yv = G_time.median_t_raw_ms(idx);   % visibility별 local time의 median

    % q 기준 정렬
    [qv, order] = sort(qv);
    yv = yv(order);

    if ~isempty(qv)
        plot(qv, yv, '-o', 'DisplayName', sprintf('\\sigma = %.3f', s));
    end
end
xlabel('visibility q');
ylabel('median t_{raw} [ms]');
%title(sprintf('Local iteration time vs visibility (iter = %d)', it_ref_time));
legend('Location','best');
hold off;
%% sdf-inlier
it_ref = 4;
T1 = T(T.iter == it_ref & T.q <= 0.6, :);

figure; grid on; hold on;
scatter(T1.inlierR, T1.SDF_RMS, 10, T1.q, 'filled');  % 색으로 q 표현
xlabel('inlier ratio');
ylabel('SDF');
ylim([0 1]);
%title('SDF\_RMS vs inlier ratio (iter = 4, q \le 0.6)');
cb = colorbar;
cb.Label.String = 'visibility q';
hold off;
%% sigma-q graph
it_ref = 4;
T1 = T(T.iter == it_ref, :);

G = groupsummary(T1, {'q','sigma'}, 'mean', 'SDF_RMS');

q_vec = sort(unique(G.q));
s_vec = sort(unique(G.sigma));

M = nan(numel(s_vec), numel(q_vec));  % [sigma, q] 매트릭스
for i = 1:numel(s_vec)
    for j = 1:numel(q_vec)
        idx = (abs(G.sigma - s_vec(i)) < 1e-12) & ...
              (abs(G.q - q_vec(j))     < 1e-12);
        if any(idx)
            M(i,j) = G.mean_SDF_RMS(idx);
        end
    end
end

figure;
imagesc(q_vec, s_vec, M);
set(gca,'YDir','normal');
xlabel('visibility q');
ylabel('\sigma (noise std)');
title('mean SDF\_RMS over (q, \sigma), iter = 4');
h = colorbar;
h.Label.String = 'mean SDF\_RMS';
%% 이전 피팅
%{
PPm_use_unbias =  PP_use - Disp;

if ~isempty(Beta)
    f1_RAN = TermsC * Beta;
    Funcs1_RAN = matlabFunction(f1_RAN);
    
    PP_inlier = P_use;
    %PP_inlier = P_use(logical(inlierMask), :);
    PP_inlier_unbias = PP_inlier - Disp;
    Value_inlier_unbias = Funcs1_RAN(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,[-20 20 -10 10 -30 30],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 150);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

    
    figure(3)
    scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,'b','filled');
    
    fimplicit3(f1_RAN-1,[-20 20 -20 20 -30 30],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor','none','FaceAlpha',0.5, 'MeshDensity', 150);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

end
%}