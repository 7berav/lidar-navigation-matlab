%2025 03 31
function [ResultDisp, ResultRot,Beta,initVal, optVal] = PoliNavigationSolver2(isGlobalApproach, PPcoord, order,binit,qinit)
    %syms x;
    %syms y;
    %syms z;

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
    
    tic
    [ResultRot,initVal, optVal] = Rotation(Beta,order,binit,qinit);
    toc
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
    tic                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    syms x y z
    TermsB = homogeneTerm(order);
    betaB = sym('beta', [1, length(TermsB)]);
    f2 = sum(betaB .* TermsB);
    FuncsB = matlabFunction(TermsB);
    toc
    tic
    difx = diff(f2, x);
    dify = diff(f2, y);
    difz = diff(f2, z);
    toc
    tic
    f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], betaB});
    d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], betaB});
    d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], betaB});
    d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], betaB});
    toc
    tic
    coord_use = coord;
    [beta_values,error]    = regressionFourthOrder( coord_use,FuncsB);
    toc

    error_shift = error;  % 초기 에러 설정
    
    shift_1= [0;0;0];
    shiftSet = [];
    shift1Set = [];
    shiftResidueSet = [0;0;0];
    shiftSum = [0 ; 0 ; 0];
    numIterations = 10;  % 반복 횟수 설정

    

    tic
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
    toc
    DispOut = shiftSum.';
    Beta = beta_values.';

end
  

function [QOut,initVal,optVal] = Rotation(Beta,order,binit,qinit)
%    syms q0 q1 q2 q3 real
%    syms x y z
    % 간단히 벡터화
%    q = [q0; q1; q2; q3];
    
    % 회전행렬 R(q) 정의
%    Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
%           2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
%           2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];
    
%    xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
%    yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
%    zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;
    

        
    %end
    
    global Mhandle
    w = ones(length(binit),1);
    w(1) = 240;
    w(2) = 240;
    w(3) = 240;
    fObj = @(qVec) sum( w .* ((Mhandle(qVec(1),qVec(2),qVec(3),qVec(4)) * Beta.' - binit).^2));
    %fNum = matlabFunction(f_coeffsB_norm, 'Vars', [q0, q1, q2, q3]);
    %fHandle3 = @(var) fNum(var(1), var(2), var(3), var(4));
    nonlcon = @(qVec) deal([],qVec'*qVec - 1);  
    

    
    global quaterseed
    %size(quaterseed);
    spset = CustomStartPointSet(quaterseed);

    %initQ = [1; 0; 0; 0];
    initQ = qinit;
    initVal = fObj(initQ.');
   

    
    
    optionsC = optimoptions('fmincon','Display','none','Algorithm','interior-point',...
        'OptimalityTolerance',3e-4,'ConstraintTolerance',1e-4,'MaxIterations',39,'UseParallel',false);
    
    problem = createOptimProblem('fmincon','objective', @(qVec) fObj(qVec), ...
    'x0', initQ, 'lb', -ones(1,4), 'ub', ones(1,4), 'nonlcon', @(qVec) nonlcon(qVec),'options', optionsC);
    
    
    
    gs = GlobalSearch('NumTrialPoints', 1650,'NumStageOnePoints',180,'StartPointsToRun',...
        'bounds','Display','none');
    gs.BasinRadiusFactor = 0.5;  % 기본값보다 작게 설정하면 basin 축소 폭이 달라짐
    gs.MaxWaitCycle = 10; 
    %gs = GlobalSearch('Display','iter','StartPointsToRun','bounds');
    %'NumTrialPoints' = 20개의 후보점을 뿌려서 유망한 지점만 local solver로 보냄
    
    [xGlobal, fvalGlobal] = run(gs, problem);
    



    %[qOpt, fValB] = fmincon(@(qIn) fObj (qIn), ...
    %                       initQ,[],[],[],[],[],[], nonlcon, optionsC);


    %qOpt = qOpt / norm(qOpt);  % safety normalize
    %optVal = fValB;

    %disp(fvalMulti);
    %disp(fvalGlobal);
    %disp(fValB);
    %disp('-----');

    qOpt = xGlobal / norm(xGlobal);
    optVal = fvalGlobal;

    QOut = qOpt;
end