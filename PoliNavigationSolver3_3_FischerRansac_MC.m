%2025 03 31
function [ResultDisp,  Beta, inlierMask,Log] = ...
         PoliNavigationSolver3_3_FischerRansac_MC(isGlobalApproach, PPcoord, order,nT,Funcs,Grads,ransacPar)
    if nargin < 5 || isempty(Funcs)
        syms x y z
        Terms = homogeneFischerTerms(order);
        Funcs  = matlabFunction(Terms, 'Vars', [x,y,z]);
        nT = numel(Terms);
    end
    
    PPm = PPcoord;                               % 원본
    center_shift = mean(PPm,1);
    PPm_shift   = PPm - center_shift;            % 원점 이동

    if isGlobalApproach
        [ResultDisp0, Beta, inlierMask] = ...
            RansacDisplacementGlobal(PPm_shift,nT,Funcs, ransacPar); %필요없음. 무시하셈
    else
        [ResultDisp0, Beta, inlierMask,Log] = ...
            RansacWeightedSingleModel(PPm_shift,nT,Funcs,Grads, ransacPar);
    end
    ResultDisp0;
    ResultDisp = ResultDisp0 + center_shift;     % 원점 복원

    % ── 이후 회전 최적화는 기존 함수 재사용 ──────────────────────────────

end

% ============================================================================
% RANSAC (Weighted Single Model)
% ============================================================================
function [bestDisp, bestBeta, bestInlierMask] = RansacDisplacementLocal(coord, nT,Funcs, p)
    % coord           : N×3,  원점 이동된 점군
    % order           : 다항식 차수
    % p               : ransacPar 구조체
    % bestDisp        : Σshift_k,   your 기존 DispOut 개념
    % bestBeta        : 최종 Beta  (cols = #term)
    % bestInlierMask  : N×1 logical

    % ── 사전 계산 ──────────────────────────────────────────────
    terms  = homogeneTerm(order);      % 아래 함수 이미 보유
    k      = 2*numel(terms);             % 최소 샘플 크기
    N      = size(coord,1);
 
    %Funcs  = matlabFunction(terms);
    % RANSAC 반복 횟수 N_i = log(1-conf) / log(1-w^k)
    if isfield(p,'minInlierRatio')
        w = p.minInlierRatio;
    else
        w = 0.5;
    end
    %p.maxIter = min(p.maxIter, ...
    %    ceil(log(1-p.conf)/log(max(realmin,1-w^k))));

    bestScore = 0;   bestBeta = [];  bestDisp = [0;0;0];
    bestInlierMask = false(N,1);
    tgt = 1;
    for iter = 1:p.maxIter
        % 1) 무작위 최소 샘플 선택
        idx = randperm(N,k);
        subset = coord(idx,:);

        % 2) 모델 추정 ─ 기존 DisplacementLocal의 “shift 최적화” 없이
        %    첫 회 귀찮으면 regressionFourthOrder로만 Beta 예비 추정:

        [beta_tmp, ~] = regressionFourthOrder(subset, Funcs);

        % 3) 전체 잔차 계산
        residuals = calcPolyResidual(coord,beta_tmp, Funcs);
        
        % 4) 인라이어 집합
        inlierMask = residuals < p.thresh;
        score      = sum(inlierMask);
        
        % 5) 최고 모델 갱신
        if score > bestScore
            bestScore       = score;
            bestInlierMask  = inlierMask;
            bestBeta        = beta_tmp;
            %fprintf('Best Score (before) : %d \t',bestScore);
            if score > w*N   
                w = score/N;  % w 업데이트 
            end   
            %disp(residuals.');
            ceil(log(1-p.conf)/log(max(realmin,1-w^k)));
            p.maxIter = min(p.maxIter, ...
                  ceil(log(1-p.conf)/log(max(realmin,1-w^k))));
            inlierCoord = coord(bestInlierMask,:);

            [bestDisp, bestBeta] = DisplacementLocal(inlierCoord, order);
            residuals_temp = calcPolyResidual(coord-bestDisp,bestBeta, Funcs);
            inlierMask_temp = residuals_temp < p.thresh;
            score      = sum(inlierMask_temp);
            inlierCoord_temp = coord(inlierMask_temp,:);
            %fprintf('(after) : %d \n',score);
            

           axis equal;

            tgt = tgt + 1;
        end
        if iter >= p.maxIter, break; end
    end

    % 6) 인라이어로 모델 재추정 (+ shift loop 활용)
    inlierCoord = coord(bestInlierMask,:);
    [bestDisp, bestBeta] = DisplacementLocal(inlierCoord, order);
    

end

function [bestDisp, bestBeta, bestInMask,Log] = RansacWeightedSingleModel(coord,nT,Funcs,Grads, p)
    % PPm_use : N×3 point cloud
    % order   : polynomial order
    % Funcs   : used terms
    % p       : struct with fields
    %            p.lambda (optional, default 0): ridge penalty for regressionFourthOrder
    %
    % Returns bestBeta, bestDisp and accumulated soft‐weight W (N×1).

    % ── Ridge lambda ──────────────────────────────────────────────
    if isfield(p, 'lambda')
        lambda = p.lambda;
    else
        lambda = 0;   % OLS (기존 동작과 동일)
    end

    % ── 초기 설정 ──────────────────────────────────────────────────
    N    = size(coord,1);
    
    k    = round(p.k * nT); % 수치 크게 할수록 데이터 많이 걸리고, 딱히 성능이 나빠지지는 않음. 낮추면 애초에 좋은 샘플을 못뽑음. 

    sigma0 = 0.15;
    sigma  = sigma0;
    if isfield(p,'minInlierRatio')
        w = p.minInlierRatio;
    else
        w = 0.5;
    end


    % ====== [MC ADD] Monte Carlo 옵션/로그 버퍼 ======
    mc.on        = isfield(p,'mc') && ~isempty(p.mc) && isfield(p.mc,'on') && p.mc.on;
    mc.var       = '';   if mc.on && isfield(p.mc,'saveVarName'), mc.var  = p.mc.saveVarName; end
    mc.file      = '';   if mc.on && isfield(p.mc,'saveMatFile'), mc.file = p.mc.saveMatFile; end
    mc.hasGT = false;
    mc.gtMask = [];
    if mc.on && isfield(p.mc,'gtMask')
        mc.gtMask = logical(p.mc.gtMask(:));
        if numel(mc.gtMask) == N
            mc.hasGT = true;
        end
    end
    mc.metric = 'none';
    if mc.on && isfield(p.mc,'metric'), mc.metric = p.mc.metric; end
    mc.keepBestPerIter = true; % 고정
    if mc.on
        %score_raw, score_loc, elapsed_ms
        IterLog = nan(p.maxIter, 6);
    end

    if ~isfield(p,'sim'), p.sim = struct(); end
    if ~isfield(p.sim,'freezeW'),     p.sim.freezeW = true;     end   % w 갱신 동결
    if ~isfield(p.sim,'freezeIter'),  p.sim.freezeIter = true;  end   % p.maxIter 동적 갱신 금지
    %if ~isfield(p.sim,'freezeSigma'), p.sim.freezeSigma = false; end   % sigma 적응 금지
    if ~isfield(p.sim,'freezeOmega'), p.sim.freezeOmega = false; end   % omega(샘플링 가중) 갱신 동결
    omega_min = 1e-6;
    
    if ~p.sim.freezeIter
            p.maxIter = min(p.maxIter, ...
        ceil(log(1-p.conf)/log(max(realmin,1-w^k)))); %신뢰구간 내 샘플링 횟수
    end

    
    % 가중치 누적용
    omega    = ones(N,1)*0.33; 
    omegajLoc= zeros(N,1);
    Ssum = 4.6;% 클수록 천천히 w 반영 /0.4 하면 대충 몇번 후에 초기가중치 탈피하는지 판단 가능. 2,3만 아니면 됨
    SjLoc = NaN;

    bestScore = 0.13; %어짜피 갱신될거 충분히 작게. 
    bestBeta  = [];
    bestDisp  = [0,0,0];
    bestInMask = false(N,1);
    target = 1;

    Log = struct( ...
        'preScore', {}, ...   % Sj
        'postScore', {}, ...  % SjLoc
        't_raw_ms', {}, ...   % raw 단계 시간[ms]
        't_loc_ms', {}, ...   % local 단계 시간[ms]
        'w', {}, ...          % 그 iter의 w
        'maxN', {}, ...       % 그 iter에서 계산된 maxIter 후보
        'inlierR', {}, ...    % inlier ratio
        'SDF_RMS', {}, ...    % 그 iter 모델 기준 SDF RMS
        'SDF_P50', {}, ...    % 그 iter 모델 기준 SDF P50
        'Disp', {}, ...
        'betaNorm', {}, ...   % ||beta||^2 — ridge 정규화 효과 확인용
        'TP', {}, ...
        'FP', {}, ...
        'FN', {}, ...
        'TN', {}, ...
        'precision', {}, ...
        'recall', {}, ...
        'F1', {} ...
    );
    for iter = 1:p.maxIter
        raw_tic = tic;
        w_for_log     = w;
        maxN_for_log  = NaN;
        local_ms   = NaN;
        SjLoc      = NaN;
        SDF_RMS    = NaN;
        SDF_P50    = NaN;
        disp_for_log = [NaN NaN NaN];
        inlierR_for_log = NaN;
        betaNorm_for_log = NaN;

        TP_for_log = NaN;
        FP_for_log = NaN;
        FN_for_log = NaN;
        TN_for_log = NaN;
        prec_for_log = NaN;
        rec_for_log  = NaN;
        F1_for_log   = NaN;

        % 1) uniform 또는 가중치 기반 샘플링
        idx = randsample(1:N, k, true, omega);%omega 기반 샘플링
        % 2) provisional β  (ridge if lambda > 0)
        coord_unbias = coord-mean(coord(idx,:),1);
        betaTmp = regressionFourthOrder(coord_unbias(idx,:), Funcs, lambda);
        r    = abs(Funcs([coord_unbias(:,1),coord_unbias(:,2),coord_unbias(:,3)])*betaTmp - 1);
        omegaj   = exp(-r.^2/(2*sigma0^2)); %이 iter의 점간 점수 0~1
        % 4) score = mean(wj)
        Sj = mean(omegaj);
        raw_ms = 1000 * toc(raw_tic);
        % 갱신 전 기록
        Log(iter).preScore = double(Sj);
        Log(iter).t_raw_ms = double(raw_ms);



        % 5) 갱신 조건: “최고 점수의 updateThresh 이상일 때”
        if Sj >= bestScore * p.updateThresh
            % 로컬 최적화 (shift + refit)
            local_tic = tic;
            inMaskj = r < p.thresh;
            %w = 0.9*w+0.1*sum(omega(inMaskj))/sum(omega);
            if ~p.sim.freezeW
                w = 0.50*w + 0.50*sum(omega(inMaskj))/sum(omega);
            end
            if ~p.sim.freezeIter
                p.maxIter = min(p.maxIter, ...
                    2*ceil(log(1-p.conf)/log(max(realmin,min(1-1e-12,1-w^k)))));
                %fprintf('%d',p.maxIter);
            end
            %Log(iter).w        = double(w);
            %Log(iter).maxN     = double(2*ceil(log(1-p.conf)/log(max(realmin,min(1-1e-12,1-w^k)))));

            [dispLoc, betaLoc] = DisplacementLocal(coord(inMaskj,:), Funcs,Grads,p,lambda);
            rLoc  = abs( Funcs([coord(:,1)-dispLoc(1),coord(:,2)-dispLoc(2),coord(:,3)-dispLoc(3)])*betaLoc - 1 );
            inMaskLoc = rLoc < p.thresh;
            wLoc = sum(omega(inMaskLoc))/sum(omega); % w
            

            omegajLoc = exp(-rLoc.^2/(2*sigma^2)); %점수 
            SjLoc = mean(omegajLoc); %점수 증분
            if ~p.sim.freezeOmega
                omega        = (Ssum / (Ssum + SjLoc))*omega + (SjLoc / (Ssum + SjLoc)) * omegajLoc;
                Ssum         = Ssum + SjLoc;    
                %fprintf('(bef): %d  \t',Sj);
                %fprintf('(aft): %d \t %d %n',SjLoc,wLoc);
            end
            if SjLoc >= bestScore

                bestScore = SjLoc;
                bestBeta  = betaLoc;
                bestDisp  = dispLoc;
                bestInMask= inMaskLoc;

                target = target+1;
            end
            disp_for_log = dispLoc(:).';
            local_ms = 1000 * toc(local_tic);
            
            %{
            P_result = [coord(:,1)-dispLoc(1), ...
                        coord(:,2)-dispLoc(2), ...
                        coord(:,3)-dispLoc(3)];      % N×3

            Phi   = Funcs(P_result);                 % N×nT
            Fval  = Phi * betaLoc - 1;              % N×1

            Gstack = Grads(P_result);               % (3N)×nT
            Nn     = size(P_result,1);
            gv     = Gstack * betaLoc;              % (3N)×1

            gn = sqrt( max(1e-12, ...
                     gv(1:Nn).^2 + gv(Nn+1:2*Nn).^2 + gv(2*Nn+1:3*Nn).^2) );

            d_sdf   = Fval ./ gn;                   % N×1
            SDF_RMS = sqrt(mean(d_sdf.^2));
            SDF_P50 = prctile(d_sdf, 50);
            %}
            %Log(iter).postScore = double(SjLoc);
            %Log(iter).t_loc_ms  = double(local_ms);
            %Log(iter).inlierR   = mean(inMaskLoc);

            inlierR_for_log     = mean(inMaskLoc);
            betaNorm_for_log    = betaLoc(:).' * betaLoc(:);   % ||beta||^2
            
            if mc.hasGT 
                gt = mc.gtMask;   % N×1 logical
        
                TP = sum( gt  &  inMaskLoc );
                FP = sum(~gt  &  inMaskLoc );
                FN = sum( gt  & ~inMaskLoc );
                TN = sum(~gt  & ~inMaskLoc );
        
                precision = TP / max(1, TP + FP);
                recall    = TP / max(1, TP + FN);
                F1        = 2 * precision * recall / max(1e-12, precision + recall);
        
                TP_for_log = TP;
                FP_for_log = FP;
                FN_for_log = FN;
                TN_for_log = TN;
                prec_for_log = precision;
                rec_for_log  = recall;
                F1_for_log   = F1;
            end
            %Log(iter).SDF_RMS   = double(SDF_RMS);
            Log(iter).SDF_P50   = double(SDF_P50);
            Log(iter).Disp      = disp_for_log;
            
           
            
        end
        if mc.on
            IterLog(iter, 1:2) = [double(Sj), double(SjLoc)];
            IterLog(iter, 3:4) = [raw_ms,local_ms];
            IterLog(iter, 5:6) = [w,p.maxIter];     
        end
        Log(iter).postScore = double(SjLoc);
        Log(iter).t_loc_ms  = double(local_ms);
        Log(iter).inlierR   = double(inlierR_for_log);
        Log(iter).betaNorm  = double(betaNorm_for_log);
        Log(iter).w        = double(w);
        Log(iter).maxN     = double(2*ceil(log(1-p.conf)/log(max(realmin,min(1-1e-12,1-w^k)))));
        Log(iter).TP        = double(TP_for_log);
        Log(iter).FP        = double(FP_for_log);
        Log(iter).FN        = double(FN_for_log);
        Log(iter).TN        = double(TN_for_log);
        Log(iter).precision = double(prec_for_log);
        Log(iter).recall    = double(rec_for_log);
        Log(iter).F1        = double(F1_for_log);

        if iter >= p.maxIter
                
                fprintf('!!break!!\n');
                break;
        end
    end
    

    %{
    inlierCoord = coord(bestInMask,:);
    figure(4);
    scatter3(coord(:,1),coord(:,2),coord(:,3),2,omega(:),'filled');
    hold on
    %scatter3(inlierCoord(:,1), inlierCoord(:,2), inlierCoord(:,3), 2, [0.5 0.5 0.5], 'filled');
    hold off
    colorbar
    axis equal;


    figure(5);
    scatter3(coord(:,1),coord(:,2),coord(:,3),2,omegajLoc(:),'filled');
    hold on
    %scatter3(inlierCoord(:,1), inlierCoord(:,2), inlierCoord(:,3), 2, 'r', 'filled');
    hold off
    colorbar
    axis equal;
    drawnow;
    %}
    if exist('mc','var') && mc.on
        % 유효 행만 남기기 (앞쪽 NaN 허용, 뒤쪽 미사용 row 제거)
        lastIter = find(~isnan(IterLog(:,1)), 1, 'last');
        if ~isempty(lastIter), IterLog = IterLog(1:lastIter, :); end
    
        % 작업공간(base)에 저장 (옵션)
        if ~isempty(mc.var)
            assignin('base', mc.var, IterLog); 
        end
    
        % MAT 파일 저장 (옵션)
        if ~isempty(mc.file)
            save(mc.file, 'IterLog');
        end
    end
end


    
function [DispOut, Beta] = DisplacementLocal(coord, Funcs, Grads, p, lambda)
    % lambda (optional): ridge penalty passed from RansacWeightedSingleModel
    if nargin < 5 || isempty(lambda)
        lambda = 0;
    end
    if nargin<4 || isempty(p)
        p.locIters  = 7;
        p.damping   = 0.85;
        p.momentum  = 0.66;
        p.reg       = 1e-6;
        p.tol       = 0;       % 0이면 미사용
        p.verbose   = false;
    end

    coord_center = mean(coord, 1);
    coord_use = coord  - coord_center;
    [beta_values,error]    = regressionFourthOrder(coord_use, Funcs, lambda);
    error_shift = error;  % 초기 에러 설정


    centerSum = [0 ; 0 ; 0];
    v        = [0;0;0];
    N  = size(coord_use,1);
    for it = 1:p.locIters
        % 1) 그래디언트: (3N×nT)*(nT×1) → (3N×1) → N×3
        Gstack = Grads(coord_use);            % (3N)×nT
        GradVal   = Gstack * beta_values(:);            % (3N)×1

        %dx = d1_numeric( coord_use, beta_values.');
        %dy = d2_numeric( coord_use, beta_values.');
        %dz = d3_numeric( coord_use, beta_values.');


        Gx = GradVal(1:N);
        Gy = GradVal(N+1:2*N);
        Gz = GradVal(2*N+1:3*N);
        Gxyz = [Gx, Gy, Gz];

        [delta, ~] = regressionShift( coord_use, error_shift,Gxyz);% 함수 결과값부터 부호가 반대


        v      = p.momentum*v + p.damping*delta;
        center  = -v ;
        centerSum = centerSum + center;



        coord_use =  coord_use - center.';

        [beta_values, error_shift] = regressionFourthOrder(coord_use, Funcs, lambda);
        
        %0.001s 
        
    end
    
    DispOut = centerSum.';
    Beta = beta_values;

end
  
function r = calcPolyResidual(coord, beta, Funcs)
% beta: column vector, length = numel(terms)
    phi = Funcs(coord(:,1),coord(:,2),coord(:,3)); % N × #term
    r   = abs(phi*beta-1);                           % N×1
end

