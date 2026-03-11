%% 공통 세팅
totalN       = 10000;
eps_list     = [0.10 0.2 0.3]; %[0.00 0.10 0.20 0.30];
k_list       = [1.4];            % nT 배수 (저장명에는 round(10*k) 사용)
order_list   = [4];
Niter        = 10000;

PP11 = generateRandomPointsOnHexagonPrism(10400) + randn(10400,3)*0.007;  % inlier pool
PP14 = (2*rand(5000,3)-1)*2;                                             % outlier pool (Uniform[-2,2])
PP12 = generateRandomPointsOnCylinder(10400) + randn(10400,3)*0.007;  % inlier pool
PP13 = generateRandomPointsOnCube(10400) + randn(10400,3)*0.007;  % inlier pool

% 고정 RANSAC 파라미터(필요 필드만 세팅)
ransacPar = struct( ...
  'maxIter', Niter, ...
  'conf', 0.96, ...
  'thresh', 0.40, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0.60, ...
  'locIters', 4, ...
  'momentum', 0.66, ...
  'damping', 0.85, ...
  'reg', 1e-6, ...
  'k', 1.4 ...              % 루프에서 덮어씀
);
k=1.4;
%ransacPar.mc  = struct('on', true, 'saveVarName', 'MC_scores_tmp', 'saveMatFile','', 'time', true, 'localOff', false);
ransacPar.mc  = struct( ...
    'on', true, ...
    'saveVarName', 'MC_scores_tmp', ...
    'saveMatFile','', ...
    'time', true, ...
    'localOff', false, ...
    'metric', 'basic' ...         % N×1 logical, ground truth inlier mask
);
ransacPar.sim = struct('freezeW', false, 'freezeIter',false, 'freezeOmega', false);

warning('off','MATLAB:nearlySingularMatrix');
%%
ResRAN = struct([]);
row    = 0;
%%
BestRAN = struct([]);
rowBest = 0; 
%%
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
nQ = 2;
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

sigma_list = [0.01];                % 노이즈 σ
q_list     = [0.30 0.0];  % 상위 q 분위수 컷(남는 비율 ≈ 1-q)
it_list    = [4];

%%
for order = order_list
    % ----- order별 basis/미분 생성 -----
    TermsC = homogeneFischerTerms(order);
    [Funcs1, Grads1, ~] = makeFuncsGradsStack(TermsC);
    nT = numel(TermsC);  
    for it = it_list
        ransacPar.locIters = it;          % 로컬 반복 수
    
        for pid = 1:5
            P0 = P_sets{pid};             % 기본 점군 (N0×3)
    
            for Qid = 1:nQ
                R = Q_list{Qid};          % 회전행렬 (3×3)
                Prot = (R * P0.').';      % 회전 적용: N0×3
    
                Pnoisy = Prot+ randn(size(Prot))*0.01;         
    
                for q = q_list
                    % ---- 가시율: z-축 q 분위수 이상만 사용 ----
                    zproj = Pnoisy(:,3);
                    tau   = quantile(zproj, q);
                    keep  = (zproj >= tau);
                    P_vis = Pnoisy(keep, :);           % 가시 점들
    
                    if size(P_vis,1) < 2*nT
                        fprintf('[skip] pid=%d, Qid=%d, q=%.2f → N=%d\n', ...
                                pid, Qid, q, size(P_vis,1));
                        continue;
                    end
    
                    % ---- 인라이어/아웃라이어 비율 eps 루프 ----
                    for eps = eps_list    % eps = outlier 비율이라고 가정
                        totalN = min(size(P_vis,1), 10000);   % 전체 사용 점 수
                        n_out  = round(totalN * eps);
                        n_in   = totalN - n_out;
    
                        if n_in <= 0
                            continue;
                        end
    
                        % inlier: P_vis 일부
                        Pin = P_vis(1:n_in, :);
    
                        % outlier: 예를 들어 [-2,2]^3 uniform
                        Pout = (2*rand(n_out,3) - 1)*2;
    
                        PPm_use = [Pin; Pout];         % N×3 (N = totalN)
    
                        % ---- GT inlier mask (길이 N) ----
                        gtMask = false(size(PPm_use,1),1);
                        gtMask(1:n_in) = true;
    
                        % ---- RANSAC 파라미터 세팅 ----
                        %ransacPar.k = k;         % 또는 따로 k_list 돌릴 거면 바깥 루프
                        ransacPar.mc.gtMask     = gtMask;
                        ransacPar.mc.metric     = 'basic';  % precision/recall 계산
                        ransacPar.mc.saveVarName = '';      % 필요 없으면 비움
    
                        totaltic = tic;
    
                        % ---- RANSAC 호출 ----
                        [DispRAN, BetaRAN, inlierMaskRAN, LogIter] = ...
                            PoliNavigationSolver3_3_FischerRansac_MC( ...
                                0, PPm_use, order, numel(TermsC), Funcs1, Grads1, ransacPar);
    
                        t_total = toc(totaltic);      % 전체 RANSAC 시간 [s]
                        disp(t_total);
                        % ---- 최종 best 모델에 대한 지표 계산 ----
                        Npts = size(PPm_use,1);
    
                        % (1) 모델 기준 inlier 비율 (SDF 기준)
                        P_result = PPm_use - DispRAN;       % 모델 좌표계
                        Phi      = Funcs1(P_result);        % N×nT
                        Fval     = Phi * BetaRAN - 1;       % N×1
                        inlierMask_model = abs(Fval) <= ransacPar.thresh;
                        inlierR_model = mean(inlierMask_model);  % 모델 기준 인라이어 비율
    
                        % (2) GT 기준 recall / precision
                        TP = sum( gtMask  &  inlierMaskRAN );
                        FP = sum(~gtMask  &  inlierMaskRAN );
                        FN = sum( gtMask  & ~inlierMaskRAN );
                        TN = sum(~gtMask  & ~inlierMaskRAN );
    
                        precision = TP / max(1, TP + FP);
                        recall    = TP / max(1, TP + FN);
                        F1        = 2 * precision * recall / max(1e-12, precision + recall);
    
                        % ---- 요약 구조체에 한 줄 기록 ----
                        rowBest = rowBest + 1;
    
                        BestRAN(rowBest).it      = it;
                        BestRAN(rowBest).pid     = pid;
                        BestRAN(rowBest).Qid     = Qid;
                        BestRAN(rowBest).q       = q;
                        BestRAN(rowBest).eps     = eps;
                        BestRAN(rowBest).order   = order;
                        BestRAN(rowBest).k       = ransacPar.k;
                        BestRAN(rowBest).N       = Npts;
    
                        BestRAN(rowBest).inlierR   = inlierR_model;
                        BestRAN(rowBest).recall    = recall;
                        BestRAN(rowBest).precision = precision;
                        BestRAN(rowBest).F1        = F1;
                        BestRAN(rowBest).t_total   = t_total;
    
                        %fprintf('it=%d pid=%d Q=%2d q=%.2f eps=%.2f N=%5d | rec=%.3f | inlierR=%.3f | t=%.3fs\n', ...
                        %    it, pid, Qid, q, eps, Npts, recall, inlierR_model, t_total);
    
                    end % eps
                end % q
            end % Qid
        end % pid
    end % it
end
T_best = struct2table(BestRAN);
%{
newOrder = { ...
    'pid', 'Qid', 'sigma', 'error', 'q', 'iter', 'N', ...   % 앞에 오게 하고 싶은 것들
    'preScore', 'postScore', 't_raw_ms', 't_loc_ms', ...
    'w', 'maxN', 'inlierR', 'SDF_RMS', 'SDF_P50', 'Disp', ...
    'TP', 'FP', 'FN', 'TN', 'precision', 'recall', 'F1' ... 
};
ResRAN = orderfields(ResRAN, newOrder);
T = struct2table(ResRAN);
%}

%%
S = load('monte_RANSAC/T_best_localafter.mat');
T = S.T_best;

%T = T_best;
figure; hold on;
for q_vis = [0.0, 0.3]
    Tq = T(abs(T.q - q_vis) < 1e-6, :); % 가시율 필터

    
    orders = unique(Tq.order);
    for oi = 1:numel(orders)
        ord = orders(oi);
        Tqo = Tq(Tq.order == ord, :);

        % eps × order 그룹별 평균/표준편차
        [G, eps_vals] = findgroups(Tqo.eps);
        F1_mean = splitapply(@mean, Tqo.F1, G);
        F1_std  = splitapply(@std,  Tqo.F1, G);

        % 단순선 그래프
        %plot(eps_vals, F1_mean, '-o', 'DisplayName', sprintf('order=%d', ord));
        h =errorbar( eps_vals, F1_mean, F1_std, '-o', 'DisplayName', sprintf('order=%d', ord));
        %set(h, 'HandleVisibility', 'off');
    end
    xlabel('Outlier ratio');
    ylabel('F1 score');
    xlim([0 0.4]);
    legend('Location','best');
    grid on;
end
%%
figure; hold on;
for q_vis = [0.0, 0.3]
    Tq = T(abs(T.q - q_vis) < 1e-6, :); % 가시율 필터

    
    orders = unique(Tq.order);
    for oi = 1:numel(orders)
        ord = orders(oi);
        Tqo = Tq(Tq.order == ord, :);

        % eps × order 그룹별 평균/표준편차
        [G, eps_vals] = findgroups(Tqo.eps);
        pre_mean = splitapply(@mean, Tqo.precision, G);
        pre_std  = splitapply(@std,  Tqo.precision, G);

        % 단순선 그래프
        %plot(eps_vals, pre_mean, '-o', 'DisplayName', sprintf('order=%d', ord));
        
         errorbar(eps_vals, pre_mean, pre_std, '-o', 'DisplayName', sprintf('order=%d', ord))
    end
    xlabel('Outlier ratio');
    ylabel('Precision');
    xlim([0 0.4]);
    legend('Location','best');
    grid on;
end

%% time
    figure; hold on;
for q_vis = [0.0, 0.3]
    Tq = T(abs(T.q - q_vis) < 1e-6, :);


    orders = unique(Tq.order);
    for oi = 1:numel(orders)
        ord = orders(oi);
        Tqo = Tq(Tq.order == ord, :);

        [G, eps_vals] = findgroups(Tqo.eps);
        t_mean = splitapply(@mean, Tqo.t_total, G);
        t_std  = splitapply(@std,  Tqo.t_total, G);

        errorbar(eps_vals, t_mean, t_std, '-o', ...
                 'DisplayName', sprintf('order=%d', ord));
    end
    xlabel('Outlier ratio');
    ylabel('Total time [s]');
    xlim([0 0.4]);
    title(sprintf('time vs outlier (q=%.2f)', q_vis));
    legend('Location','best');
    grid on;
end

%%
%{
%% 시계열분석


target_eps   = 0.20;
%target_k     = 1.2;
target_order = 4;

idx = T.error == target_eps ;

Tsub = T(idx,:);
iterIdx = (1:height(Tsub)).';
figure('Position',[100 100 560 740]);
tiledlayout(3,1);

% 1) w
nexttile;
plot(iterIdx, Tsub.w, '-','LineWidth', 1.2);
xlabel('Iteration');
ylabel('w');
%title(sprintf('w vs iter (eps=%.2f, k=%.2f, order=%d)', target_eps, target_k, target_order));
grid on;

% 2) maxIter (maxN)
nexttile;
plot(iterIdx, Tsub.maxN, '-','LineWidth', 1.2);
xlabel('Iteration');
ylabel('maxN (estimated maxIter)');
xlim([0 10000]); 
ylim([0 20000]);
%title('maxIter estimate vs iter');
grid on;

% 3) inlier ratio
nexttile;
plot(iterIdx, Tsub.inlierR, '.-','MarkerSize', 6,'LineWidth', 1.2);
xlabel('Iteration');
ylabel('inlierR');
ylabel('maxN (estimated maxIter)');
xlim([0 10000]); 
ylim([0 1])
%title('inlier ratio vs iter');
grid on;

%% 정밀도
figure('Position',[100 100 560 740]);
tiledlayout(3,1);

% 1) precision
nexttile;
plot(iterIdx, Tsub.precision, '.-', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('Iteration'); ylabel('Precision');
xlim([0 10000]);
ylim([0 1]); grid on;

% 2) recall
nexttile;
plot(iterIdx, Tsub.recall, '.-', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('Iteration'); ylabel('Recall');
xlim([0 10000]);
ylim([0 1]); grid on;

% 3) F1
nexttile;
plot(iterIdx, Tsub.F1, '.-', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('Iteration'); ylabel('F1');
xlim([0 10000]);
ylim([0 1]); grid on;

%}
