%2024 08 20
tic
PP01 =generateRandomPointsOnSurface(2403)+randn(2403,3)*0.0001;
PP02 =generateRandomPointsOnCylinder(170)+randn(170,3)*0.03;
toc
shiftReal =  [.10 -0.27 -0.20];
PPm = PP02 + shiftReal;
PPmm= PP02*3 + shiftReal;


PP2=PP02(PP02(:,3)>0.5,:);
PP3=PP02(PP02(:,1)>0.80|PP02(:,1)<-.70,:);
PP4=PP02(PP02(:,2)>0.72&PP02(:,3)>0,:);
PP5=PP02(PP02(:,2)<-0.73&PP02(:,3)<0,:);
PP6=PP02(PP02(:,2)>0.652&PP02(:,3)<0&PP02(:,1)>0,:);

PPm2= [PP2] + shiftReal;
PPm24= [PP2;PP4] + shiftReal;
PPm3= [PP2;PP4]*3 + shiftReal;
PPm4= [PP4] + shiftReal;

%R_PP2 = PP * rotationM';
 
%{
figure();
hold on
scatter3(PP2(:,1),PP2(:,2),PP2(:,3),'red');
scatter3(PP3(:,1),PP3(:,2),PP3(:,3),'b');
scatter3(PP4(:,1),PP4(:,2),PP4(:,3),'green');
scatter3(PP5(:,1),PP5(:,2),PP5(:,3),'black');
scatter3(PP6(:,1),PP6(:,2),PP6(:,3),'y');
hold off
%}


%%
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


[beta_value0,error]    = regressionFourthOrder(PP02);
PPm_use = PPm;
[beta_values,error]    = regressionFourthOrder(PPm_use);
%beta 고정
%0.00026s
tic
f2_partial = @(xyz) f2_numeric(xyz, beta_values.');
f2_partialgraph = @(x,y,z) f2_numeric([x,y,z], beta_values.');
toc
%개수따라 다르지만 0.0028~0.005s
%tic
%dx_test = d1_partial(PP);
%toc

%0.003
tic 
dx_test = d1_numeric(PPm_use, beta_values.');
dy_test = d2_numeric(PPm_use, beta_values.');
dz_test = d3_numeric(PPm_use, beta_values.');
toc


% 초기 값 설정
PPm_shift = PPm_use;  % 초기 PPm 설정
error_shift = error;  % 초기 에러 설정
shiftSet = [];
shiftResidue = shiftReal.';
shiftResidueSet = [];
numIterations = 150;  % 반복 횟수 설정

figure
hold on
fimplicit3(f2_partialgraph,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),'r');
hold off
axis equal

%%
for i = 1:numIterations
    % 선형 회귀 및 4차식 피팅
    dx = d1_numeric(PPm_shift, beta_values.');
    dy = d2_numeric(PPm_shift, beta_values.');
    dz = d3_numeric(PPm_shift, beta_values.');
    dxyz = [dx dy dz];
    [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);

    %소요시간 0.0008s


    
    shiftSet = [shiftSet shift];
    shiftResidue = shiftResidue + shift;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    PPm_shift = PPm_shift + shift.';

    [beta_values, error_shift] = regressionFourthOrder(PPm_shift);

    %0.001s 
    
end

f2_partialgraph = @(x,y,z) f2_numeric([x,y,z], beta_values.');
figure
hold on
fimplicit3(f2_partialgraph,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),'r');
hold off
axis equal

plot(shiftResidueSet(1,:))
figure
plot(shiftResidueSet(2,:))