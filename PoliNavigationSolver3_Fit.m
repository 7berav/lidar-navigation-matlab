function [ResultDisp, Beta, inlierMask, Log] = ...
    PoliNavigationSolver3_Fit(PPcoord, order, nT, Funcs, Grads, p)
% 회귀 1회 + DisplacementLocal 정련을 샘플만 바꿔가며 MC 수행
    % ---- 기본값 처리 ----
    if nargin < 4 || isempty(Funcs)
        syms x y z
        Terms = homogeneFischerTerms(order);
        Funcs  = matlabFunction(Terms, 'Vars', [x,y,z]);
        nT = numel(Terms);
    end

    if nargin < 6|| isempty(p), p = struct(); end
    if ~isfield(p,'scoreType'),  p.scoreType  = 'similarity'; end
    if ~isfield(p,'locIters'),   p.locIters   = 4; end
    if ~isfield(p,'damping'),    p.damping    = 0.85; end
    if ~isfield(p,'momentum'),   p.momentum   = 0.66; end
    if ~isfield(p,'reg'),        p.reg        = 1e-6; end
    if ~isfield(p,'sigma'),  p.sigma  = 0.13; end
    if ~isfield(p,'thresh'), p.thresh = 0.400; end

    assert(isfield(p,'thresh') && ~isempty(p.thresh), ...
        'fitPar.thresh를 설정하세요.');
    
    PPm = PPcoord;                               % 원본
    center_shift = mean(PPm,1);
    PPm_shift   = PPm - center_shift;   
    N = size(PPcoord,1);
    t_raw = tic;
    beta0 = regressionFourthOrder(PPm_shift, Funcs);
    
    r0    = abs(Funcs([PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3)])*beta0 - 1);     
    if strcmpi(p.scoreType,'similarity')
        score0 = mean( exp(-(r0.^2)/(2*p.sigma^2)) );
    else
        score0 = mean(r0 < p.thresh);
    end
    raw_ms = 1000*toc(t_raw);
    
    t_loc = tic;
    [dispLoc, betaLoc] = DisplacementLocal(PPm_shift, Funcs, Grads, p);
    r     = abs(Funcs([PPm_shift(:,1)-dispLoc(1),PPm_shift(:,2)-dispLoc(2),PPm_shift(:,3)-dispLoc(3)])*betaLoc - 1 );
            
    if strcmpi(p.scoreType,'similarity')
        w = exp(-(r.^2)/(2*p.sigma^2));
        score = mean(w);
        inlierMask = (r < p.thresh);
    else
        inlierMask = (r < p.thresh);
        score = mean(inlierMask);
    end
    local_ms = 1000*toc(t_loc);
    % 5) 결과(원점 복원)
    ResultDisp = dispLoc + center_shift;   % 원점 복원
    Beta       = betaLoc(:);
    Log        = [score0, score, raw_ms, local_ms];   % 1×4 벡터 (preScore, postScore, t_raw(ms), t_local(ms))
end
    
function [DispOut, Beta] = DisplacementLocal(coord, Funcs, Grads, p) %DisplacementLocal(coord,order,Funcs)
                   
    if nargin<4 || isempty(p)
        p.locIters  = 4;
        p.damping   = 0.85;
        p.momentum  = 0.66;
        p.reg       = 1e-6;
        p.tol       = 0;       % 0이면 미사용
        p.verbose   = false;
    end
    %syms x y z
    %TermsB = homogeneFischerTerms(order);
    %betaB = sym('beta', [1, length(TermsB)]);
    %f2 = sum(betaB .* TermsB); % 4변수 심볼릭
    %FuncsB = matlabFunction(TermsB); % 함수화

    %difx = diff(f2, x); % 3변수 미분
    %dify = diff(f2, y);
    %difz = diff(f2, z);
    %f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], betaB});
    %d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], betaB});
    %d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], betaB});
    %d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], betaB});
    

    coord_use = coord;
    [beta_values,error]    = regressionFourthOrder( coord_use,Funcs);
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
        
        [beta_values, error_shift] = regressionFourthOrder(coord_use,Funcs);
        
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
