%2025 03 31
function [ResultDisp, ResultRot,Beta] = PoliNavigationSolver(isGlobalApproach, PPcoord, order,qinit)
    syms x;
    syms y;
    syms z;

    PPm_use = PPcoord;

    center_shift = mean(PPm_use,1);
    PPm_shift = PPm_use - center_shift;

    if isGlobalApproach
        [ResultDisp0, Beta] = DisplacementGlobal(PPm_shift,order);
        ResultDisp = ResultDisp0 + center_shift;
    else
        [ResultDisp0, Beta] = DisplacementLocal(PPm_shift,order);
        ResultDisp = ResultDisp0 + center_shift;
    end

    ResultRot = Rotation(Beta,order,qinit);

end

function [DispOut, Beta] =  DisplacementGlobal(coord,order)
    TermsA = nonhomogeneTerm(order);
    betaA = sym('beta', [1, length(TermsA)]);
    f1 = sum(betaA .* TermsA);
    FuncsA = matlabFunction(TermsA);

    coord_use = coord; 
    [beta_values0,error0]    = regressionFourthOrder(coord_use,FuncsA);
    ErrorValue0= norm(error0,1);
    f1_total   = subs(f1,betaA,beta_values0.');
    syms x y z
    syms a b c 
    f1_shift = subs(f1_total, [x, y, z], [x+a, y+b, z+c]);
    [coeffsA, monomialA] = coeffs(f1_shift , [x,y,z]);
    
    coeffsA_unused = sym([]);
    for k = 1:length(monomialA)
        deg = feval(symengine, 'degree', monomialA(k), x)+feval(symengine, 'degree', monomialA(k), y)+feval(symengine, 'degree', monomialA(k), z);
        if deg ~= order & deg ~= 0
           coeffsA_unused(end+1) = coeffsA(k);  
        end
    end
    f_coeffsA_norm = sum(coeffsA_unused.^2);  % f_obj(a,b,c)
    
    %최적화
    initGuess = [0,0,0];
    gradF = gradient(f_coeffsA_norm, [a, b, c]);  
    fNum = matlabFunction(gradF, 'Vars', [a, b, c]);
    fHandle2 = @(var) fNum(var(1), var(2), var(3));
    optionsB = optimoptions('fsolve', ...
        'Display', 'none', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 3000);

    [xSol_A, fval_A, exitflag_B, output_B] = fsolve(fHandle2, initGuess, optionsB);
    minVal = double(subs(f_coeffsA_norm, [a,b,c], xSol_A));
    coord_shiftA = coord_use - xSol_A;

    [DispOut, Beta] = DisplacementLocal(coord_shiftA,order);
    DispOut = DispOut + xSol_A;

end
    
function [DispOut, Beta] = DisplacementLocal(coord,order)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    syms x y z
    TermsB = homogeneTerm(order);
    betaB = sym('beta', [1, length(TermsB)]);
    f2 = sum(betaB .* TermsB);
    FuncsB = matlabFunction(TermsB);

    difx = diff(f2, x);
    dify = diff(f2, y);
    difz = diff(f2, z);
    f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], betaB});
    d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], betaB});
    d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], betaB});
    d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], betaB});
    

    coord_use = coord;
    [beta_values,error]    = regressionFourthOrder( coord_use,FuncsB);
    
    error_shift = error;  % 초기 에러 설정
    
    shift_1= [0;0;0];
    shiftSet = [];
    shift1Set = [];
    shiftResidueSet = [0;0;0];
    shiftSum = [0 ; 0 ; 0];
    numIterations = 5;  % 반복 횟수 설정

    

   
    for i = 1:numIterations
        % 미분함수에 정의
        
        dx = d1_numeric( coord_use, beta_values.');
        dy = d2_numeric( coord_use, beta_values.');
        dz = d3_numeric( coord_use, beta_values.');
        
        dxyz = [dx dy dz];
        
        [shift, residual] = regressionShift( coord_use, error_shift,dxyz);% 함수 결과값부터 부호가 반대
        
        
        shift_1 = shift*0.85 + shift_1*0.66;
        shiftSum = shiftSum - shift_1;
        
        
        
        coord_use =  coord_use + shift_1.';
        
        [beta_values, error_shift] = regressionFourthOrder(coord_use,FuncsB);
        
        %0.001s 
        
    end
    
    DispOut = shiftSum.';
    Beta = beta_values.';

end
  

function QOut = Rotation(Beta, order,qinit)
    syms q0 q1 q2 q3 real
    syms x y z
    % 간단히 벡터화
    q = [q0; q1; q2; q3];
    
    % 회전행렬 R(q) 정의
    Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
           2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
           2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];
    
    xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
    yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
    zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;
    
    TermsB = homogeneTerm(order);
    f2 = Beta *  TermsB(:);
    f4_rotated = subs(f2, [x,y,z], [xr,yr,zr]);%중심에 있는 f3를 회전
    
    [coeffsB, monomialB] = coeffs(f4_rotated , [x,y,z]);
    coeffsB_unused = sym([]);
    monList = [
      0 6 0;
      2 4 0;   
      4 2 0;
      6 0 0;   
      0 0 6;
    ];
    
    for k = 1:length(monomialB) %안쓰는 항들의 계수 norm 구하기
        degx = feval(symengine, 'degree', monomialB(k), x);
        degy = feval(symengine, 'degree', monomialB(k), y);
        degz = feval(symengine, 'degree', monomialB(k), z);
        if ~ismember([degx, degy, degz], monList, 'rows');
           coeffsB_unused(end+1) = coeffsB(k);  
        end
    end
    
    f_coeffsB_norm = sum(coeffsB_unused.^2);  % f_obj(a,b,c)
    
    fNum = matlabFunction(f_coeffsB_norm, 'Vars', [q0, q1, q2, q3]);
    fHandle3 = @(var) fNum(var(1), var(2), var(3), var(4));
    nonlcon = @(qVec) deal([],qVec'*qVec - 1);  
    
    
    %initQ = [1; 0; 0; 0];
    initQ = qinit;
    optionsC = optimoptions('fmincon','Display','none','Algorithm','interior-point',...
        'OptimalityTolerance',3e-4,'ConstraintTolerance',1e-4,'MaxIterations',12);
    
    [qOpt, fValB] = fmincon(@(qIn) fHandle3 (qIn), ...
                           initQ,[],[],[],[],[],[], nonlcon, optionsC);
    qOpt = qOpt / norm(qOpt);  % safety normalize
        
    %대입
    f4_substituted = subs(f4_rotated, [q0,q1,q2,q3], qOpt.');%중심에 있는 f3를 회전
    [coeffsB2, monomialB2] = coeffs(f4_substituted , [x,y,z]);
    coeffsB2 = double(coeffsB2);
    
    QOut = qOpt;
end