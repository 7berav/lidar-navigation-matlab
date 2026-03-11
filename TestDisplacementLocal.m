    
function [DispOut, Beta] = TestDisplacementLocal(coord, Funcs, Grads, p) %DisplacementLocal(coord,order,Funcs)
                   
    if nargin<4 || isempty(p)
        p.locIters  = 12;
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
    [beta_values,error_init]    = regressionFourthOrder( coord_use,Funcs);
    error_shift = error_init;  % 초기 에러 설정
    
    
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
    
    DispOut = centerSum;
    Beta = beta_values;

end