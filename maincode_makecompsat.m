%2024 08 20
PP01 =generateRandomPointsOnHexagonPrism(2000)+randn(2000,3)*0.0005;
PP02 =generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 =generateRandomPointsOnCylinder(200)+randn(200,3)*0.0005;

shiftReal =  [0.20 -0.1 -0.1];
PPm_body = PP01 ;
PPm_body(:,3) = PPm_body(:,3) * 1.2;
PPm_panel1 = PP02 ;
PPm_panel1(:,1) = PPm_panel1(:,1) / 16;
PPm_panel1(:,2) = PPm_panel1(:,2) * 1.63;
PPm_panel1(:,3) = PPm_panel1(:,3) / 1.3;
PPm_panel1 = PPm_panel1 + [0 3.5 0];
PPm_panel2 = PP02 ;
PPm_panel2(:,1) = PPm_panel2(:,1) / 16;
PPm_panel2(:,2) = PPm_panel2(:,2) * 1.63;
PPm_panel2(:,3) = PPm_panel2(:,3) / 1.3;
PPm_panel2 = PPm_panel2 - [0 3.5 0];
%PPmm = PPmm + shiftReal;
PPm_total=[PPm_body ; PPm_panel1; PPm_panel2];


q= [cos(25/57.92) 0 sin(25/57.92) 0];
rotm= quat2rotm(q);
PPmR_body = PPm_body * rotm.';
figure
hold on
scatter3(PPm_body(:,1),PPm_body(:,2),PPm_body(:,3),'r');
scatter3(PPm_panel1(:,1),PPm_panel1(:,2),PPm_panel1(:,3),'b');
scatter3(PPm_panel2(:,1),PPm_panel2(:,2),PPm_panel2(:,3),'b');
hold off
axis equal
title('PC Model')
figure
hold on
scatter3(PPmR_body(:,1),PPmR_body(:,2),PPmR_body(:,3),'r');

hold off
axis equal
title('PC Model')
%%
syms x;
syms y;
syms z;
order = 6;
TermsA = nonhomogeneTerm(order);
betaA = sym('beta', [1, length(TermsA)]);
f1 = sum(betaA .* TermsA);

PPm_use = PPmR_body;
[beta_values0,error0]    = regressionFourthOrder(PPm_use,TermsA);
f1_numeric = matlabFunction(f1, 'Vars', {[x, y, z], betaA});
f1_partial = @(x,y,z) f1_numeric([x,y,z], beta_values0.')-1;
f1_total   = f1()
PPm_shift = PPm_use;  % 초기 PPm 설정
shiftResidue = shiftReal.';

%shift 결정하세요
%그리고 빼고 Residue 랑 PPm 조정하세요
%{
figure
hold on
fimplicit3(f1_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
hold off
axis equal
title('Initial')
%}
figure

hold on
fimplicit3(f1_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),'r');
hold off
axis equal
title('Initial')

shifted_f = subs(poly, [x, y, z], [x - a, y - b, z - c])
% 초기 값 설정
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

TermsB = homogeneTerm(order);
betaB = sym('beta', [1, length(TermsB)]);
f2 = sum(betaB .* TermsB);

difx = diff(f2, x);
dify = diff(f2, y);
difz = diff(f2, z);
f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], betaB});
d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], betaB});
d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], betaB});
d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], betaB});

[beta_values,error]    = regressionFourthOrder(PPm_use,TermsB);

error_shift = error;  % 초기 에러 설정

shiftResidueSet = [];
numIterations = 40;  % 반복 횟수 설정



f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
%%f2_partialgraph = @(x,y,z) f2_numeric([x,y,z], beta_values.');



figure
hold on
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),'r');
hold off
axis equal
title('Before')

%% 
% 회전탐지





%%
for i = 1:numIterations
    % 미분함수에 정의
    dx = d1_numeric(PPm_shift, beta_values.');
    dy = d2_numeric(PPm_shift, beta_values.');
    dz = d3_numeric(PPm_shift, beta_values.');
    dxyz = [dx dy dz];

    [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);
    shiftResidue = shiftResidue + shift;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    
    PPm_shift = PPm_shift + shift.';
    [beta_values, error_shift] = regressionFourthOrder(PPm_shift,TermsB);

    %0.001s 
    
end

f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
figure
hold on
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),'black');
%scatter3(PP02(:,1),PP02(:,2),PP02(:,3),'b');
hold off
axis equal
xlim([-4 4]);
ylim([-4 4]);
title('After')


norm(error_shift,1)
evaluateModel(PPm_shift,f2_numeric, beta_values)

figure
subplot(1,3,1)
plot(shiftResidueSet(1,:))
subplot(1,3,2)
plot(shiftResidueSet(2,:))
subplot(1,3,3)
plot(shiftResidueSet(3,:))