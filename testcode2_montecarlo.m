%2024 08 21
SimulationShiftSet = [];

syms x;
syms y;
syms z;

beta = sym('beta', [1, 15]);  % 15개의 심볼릭 변수로 beta 정의
f2= beta(1)*x.^4 + beta(2)*y.^4 + beta(3)*z.^4 + ...
    beta(4)*x.^2.*y.^2 + beta(5)*x.^2.*z.^2 + beta(6)*y.^2.*z.^2 + ...
    beta(7)*(x.^3).*y + beta(8)*(x.^3).*z + beta(9)*(y.^3).*x + ...
    beta(10)*(y.^3).*z + beta(11)*(z.^3).*x + beta(12)*(z.^3).*y - 1 + ...
    beta(13)*(x.^2).*y.*z + beta(14)*(y.^2).*z.*x + beta(15)*(z.^2).*x.*y;
%함수 세팅을 1 뺀꼴로 하면 힘들것 같다. 나중에 형상 피팅할떄 힘들다. 지금은 수동이지만. 
difx = diff(f2, x);
dify = diff(f2, y);
difz = diff(f2, z);

%0.08s
f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], beta});
d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], beta});
d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], beta});
d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], beta});

tic
for i = 1:15000
    %tic
    PP0 =generateRandomPointsOnCylinder(403)+randn(403,3)*0.01;
    
    shiftReal =  [0.11 -0.04 -0.09];
    PPm = PP0 + shiftReal;
    PP2=PP0(PP0(:,3)>0.5,:);
    PPm2= [PP2] + shiftReal;
 
    
    PPm_shift = PPm;  % 초기 PPm 설정
    [beta_values,error]    = regressionFourthOrder(PPm_shift);

    error_shift = error;  % 초기 에러 설정
    shiftResidue = shiftReal.';
    numIterations = 60;  % 반복 횟수 설정

    for i = 1:numIterations
    % 선형 회귀 및 4차식 피팅
        dx = d1_numeric(PPm_shift, beta_values.');
        dy = d2_numeric(PPm_shift, beta_values.');
        dz = d3_numeric(PPm_shift, beta_values.');
        dxyz = [dx dy dz];
        [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);

        shiftResidue = shiftResidue + shift;
        PPm_shift = PPm_shift + shift.';
        [beta_values, error_shift] = regressionFourthOrder(PPm_shift);
        
        % 익명 함수 정의

    end
    SimulationShiftSet = [SimulationShiftSet shiftResidue];
    %toc
end
toc
figure
histogram(SimulationShiftSet(1,:));
figure
histogram(SimulationShiftSet(2,:));
figure
histogram(SimulationShiftSet(3,:));