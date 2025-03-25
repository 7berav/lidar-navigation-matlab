%2025 1 14
PP02 =generateRandomPointsOnCylinder(6000)+randn(6000,3)*0.01;

shiftReal =  [1.20 -3.5 -3.1];
PPm = PP02 + shiftReal;
PPmm= PP02;
PPmm(:,3) = PPmm(:,3) * 2;
PPmm = PPmm + shiftReal;


PP2=PP02(PP02(:,2)>0.5,:);
PP6=PP02(PP02(:,3)<0&PP02(:,1)>0,:);
PP7=PP02(PP02(:,2)*1.7+PP02(:,3)>-0.2,:);
PPm2= [PP2] + shiftReal;
PPm26= [PP2;PP6] + shiftReal;
PPm7= [PP7] + shiftReal;


%% 비동차항 먼저 추정
syms x;
syms y;
syms z;
order = 4;
TermsA = nonhomogeneTerm(order);
betaA = sym('beta', [1, length(TermsA)]);
f1 = sum(betaA .* TermsA);

PPm_use = PPm7;
[beta_values0,error0]    = regressionFourthOrder(PPm_use,TermsA);
f1_numeric = matlabFunction(f1, 'Vars', {[x, y, z], betaA});
f1_partial = @(x,y,z) f1_numeric([x,y,z], beta_values0.')-1;

PPm_shift = PPm_use;  % 초기 PPm 설정

%{
figure
hold on
fimplicit3(f1_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),2,'r','filled');
hold off
axis equal
title('Initial')
%}
%% 비동차항을 동차항으로 추정
%f1_total   = subs(f1,betaA,beta_values0);
%syms a b c 
%f1_t_shift = subs(f1_total, [x, y, z], [x-a, y-b, z-c])


%% 그냥 좌표 평균 찾기 겸 점 초기세팅

center_shift = median(PPm_use,1);
PPm_shift = PPm_use - center_shift;

shiftResidue = shiftReal.'
shiftResidue = shiftResidue - center_shift.'
%%
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

[beta_values,error]    = regressionFourthOrder(PPm_shift,TermsB);

error_shift = error;  % 초기 에러 설정

shiftSet = [];
shift1Set = [];
shiftResidueSet = [];
numIterations = 40;  % 반복 횟수 설정

f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
%{
figure
hold on
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),2,'k','filled');
hold off
axis equal
view([1,1,1]);
%}

%%
shift_0= [];
shift_1= [0;0;0];
for i = 1:numIterations
    % 미분함수에 정의
    dx = d1_numeric(PPm_shift, beta_values.');
    dy = d2_numeric(PPm_shift, beta_values.');
    dz = d3_numeric(PPm_shift, beta_values.');
    dxyz = [dx dy dz];

    [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);
    shift_1 = shift*0.85 + shift_1*0.5;
    shiftSet = [shiftSet shift];
    shift1Set = [shift1Set shift_1];
    shiftResidue = shiftResidue + shift_1;% 부호가 반대인데 알아서 해석하세요
    shiftResidueSet = [shiftResidueSet shiftResidue];
    
    PPm_shift = PPm_shift + shift_1.'; %따라서 여기도 부호가 반대임
    [beta_values, error_shift] = regressionFourthOrder(PPm_shift,TermsB);
    
    %0.001s 
    
end

f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
figure
hold on
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),2,'k','filled');
%scatter3(PP02(:,1),PP02(:,2),PP02(:,3),'b');
hold off
axis equal
view([1,1,1]);

title('After')


%norm(error_shift,1)
%evaluateModel(PPm_shift,f2_numeric, beta_values)

figure
subplot(1,3,1)
plot(shiftResidueSet(1,:))
subplot(1,3,2)
plot(shiftResidueSet(2,:))
subplot(1,3,3)
plot(shiftResidueSet(3,:))

%shift 추이
figure
subplot(1,3,1)
hold on
plot(shiftSet(1,:))
plot(shift1Set(1,:))
hold off
subplot(1,3,2)
hold on
plot(shiftSet(2,:))
plot(shift1Set(2,:))
hold off
subplot(1,3,3)
hold on
plot(shiftSet(3,:))
plot(shift1Set(3,:))
hold off
