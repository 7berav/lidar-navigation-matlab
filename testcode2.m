%2024 08 20
tic
PP =generateRandomPointsOnSurface(2403)+randn(2403,3)*0.0001;
toc
shiftReal =  [.00 -0.07 -0.00];
PPm = PP + shiftReal;
PPmm= PP*3 + shiftReal;


PP2=PP(PP(:,3)>0.5,:);
PP3=PP(PP(:,1)>0.80|PP(:,1)<-.70,:);
PP4=PP(PP(:,2)>0.72&PP(:,3)>0,:);
PP5=PP(PP(:,2)<-0.73&PP(:,3)<0,:);
PP6=PP(PP(:,2)>0.652&PP(:,3)<0&PP(:,1)>0,:);

PPm2= [PP2] + shiftReal;
PPm3= [PP2;PP4] + shiftReal;
PPmm3= PPm3*3 + shiftReal;
PPm4= [PP4] + shiftReal;
q= [cos(10/180*pi) 0 sin(10/180*pi)*1/sqrt(5) sin(10/180*pi)*2/sqrt(5)];
rotationM = quat2rotm(q);
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

%% 일반 회귀 후 plot
%{
%[beta0 error0]=regressionFourthOrder([PPm])
[beta1 error1]=regressionFourthOrder([PP2;PP4]);
[beta2 error2]=regressionFourthOrder([PP4]);
[beta3 error3]=regressionFourthOrder([PP]);
[beta4 error4]=regressionFourthOrder([R_PP2]);

f=@(x,y,z) beta2(1)*x.^4+  beta2(2)*y.^4+    beta2(3)*z.^4+    beta2(4)*x.^2*y.^2   +beta2(5)*x.^2*z.^2   +beta2(6)*y.^2*z.^2 ...
   +beta2(7)*(x.^3).*y    +beta2(8)*(x.^3).*z    +beta2(9)*(y.^3).*x    +beta2(10)*(y.^3).*z    +beta2(11)*(z.^3).*x    +beta2(12)*(z.^3).*y -1 ...
   +beta2(13)*(x.^2).*y.*z+beta2(14)*(y.^2).*z.*x+beta2(15)*(z.^2).*x.*y;
g=@(x,y,z) beta1(1)*x.^4+  beta1(2)*y.^4+    beta1(3)*z.^4+    beta1(4)*x.^2*y.^2   +beta1(5)*x.^2*z.^2   +beta1(6)*y.^2*z.^2 ...
   +beta1(7)*(x.^3).*y    +beta1(8)*(x.^3).*z    +beta1(9)*(y.^3).*x    +beta1(10)*(y.^3).*z    +beta1(11)*(z.^3).*x    +beta1(12)*(z.^3).*y -1 ...
   +beta1(13)*(x.^2).*y.*z+beta1(14)*(y.^2).*z.*x+beta1(15)*(z.^2).*x.*y;

%h=@(x,y,z) Ans3(1)*x.^4+  Ans3(2)*y.^4+    Ans3(3)*z.^4+    Ans3(4)*x.^2*y.^2   +Ans3(5)*x.^2*z.^2   +Ans3(6)*y.^2*z.^2 ...
%   +Ans3(7)*(x.^3).*y    +Ans3(8)*(x.^3).*z    +Ans3(9)*(y.^3).*x    +Ans3(10)*(y.^3).*z    +Ans3(11)*(z.^3).*x    +Ans3(12)*(z.^3).*y -1 ;


hr=@(x,y,z) beta4(1)*x.^4+ beta4(2)*y.^4+    beta4(3)*z.^4+    beta4(4)*x.^2*y.^2   +beta4(5)*x.^2*z.^2   +beta4(6)*y.^2*z.^2 ...
   +beta4(7)*(x.^3).*y    +beta4(8)*(x.^3).*z    +beta4(9)*(y.^3).*x    +beta4(10)*(y.^3).*z    +beta3(11)*(z.^3).*x    +beta4(12)*(z.^3).*y -1 ...
   +beta4(13)*(x.^2).*y.*z+beta4(14)*(y.^2).*z.*x+beta4(15)*(z.^2).*x.*y;
h=@(x,y,z) beta3(1)*x.^4+  beta3(2)*y.^4+    beta3(3)*z.^4+    beta3(4)*x.^2*y.^2   +beta3(5)*x.^2*z.^2   +beta3(6)*y.^2*z.^2 ...
   +beta3(7)*(x.^3).*y    +beta3(8)*(x.^3).*z    +beta3(9)*(y.^3).*x    +beta3(10)*(y.^3).*z    +beta3(11)*(z.^3).*x    +beta3(12)*(z.^3).*y -1 ...
   +beta3(13)*(x.^2).*y.*z+beta3(14)*(y.^2).*z.*x+beta3(15)*(z.^2).*x.*y;
%}
%%
%{
beta2
figure;
hold on
fimplicit3(g,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP2(:,1),PP2(:,2),PP2(:,3),'b');
scatter3(PP4(:,1),PP4(:,2),PP4(:,3),'g');
hold off
%}


%%
%{
figure;
hold on
%fimplicit3(h,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
fimplicit3(hr,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP(:,1),PP(:,2),PP(:,3),'g');
scatter3(R_PP2(:,1),R_PP2(:,2),R_PP2(:,3),'b');

hold off

figure;
hold on
%fimplicit3(h,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
fimplicit3(h_2,[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
scatter3(PP(:,1),PP(:,2),PP(:,3),'g');


hold off
%}

%%
PPm_use = PPmm3;
[beta,error]    = regressionFourthOrder(PPm_use);
betaError = norm(beta(4:15),1)
% 초기 값 설정
PPm_shift = PPm_use;  % 초기 PPm 설정
error_shift = error;  % 초기 에러 설정

colors = {'r', 'g', 'black'};  % 색상 설정

shiftSet = [];
shiftResidue = shiftReal.';
shiftResidueSet = [];
numIterations = 150;  % 반복 횟수 설정

syms x;
syms y;
syms z;
f2= beta(1)*x.^4 + beta(2)*y.^4 + beta(3)*z.^4 + ...
    beta(4)*x.^2.*y.^2 + beta(5)*x.^2.*z.^2 + beta(6)*y.^2.*z.^2 + ...
    beta(7)*(x.^3).*y + beta(8)*(x.^3).*z + beta(9)*(y.^3).*x + ...
    beta(10)*(y.^3).*z + beta(11)*(z.^3).*x + beta(12)*(z.^3).*y - 1 + ...
    beta(13)*(x.^2).*y.*z + beta(14)*(y.^2).*z.*x + beta(15)*(z.^2).*x.*y;
%f = @(x, y, z) beta(1)*x.^4 + beta(2)*y.^4 + beta(3)*z.^4 + ...
%               beta(4)*x.^2.*y.^2 + beta(5)*x.^2.*z.^2 + beta(6)*y.^2.*z.^2 + ...
%               beta(7)*(x.^3).*y + beta(8)*(x.^3).*z + beta(9)*(y.^3).*x + ...
%               beta(10)*(y.^3).*z + beta(11)*(z.^3).*x + beta(12)*(z.^3).*y - 1 + ...
%               beta(13)*(x.^2).*y.*z + beta(14)*(y.^2).*z.*x + beta(15)*(z.^2).*x.*y;
figure
hold on
f2_numeric = matlabFunction(f2,'Vars',[x,y,z]);
fimplicit3(f2_numeric,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),'b');

hold off
axis equal

%test


%testg


%branching test
%%
for i = 1:numIterations
    % 선형 회귀 및 4차식 피팅
    1
    
    [shift, residual] = regressionShift(PPm_shift, error_shift,f2);
    
    %소요시간 0.12초 내외
    
    shiftSet = [shiftSet shift];
    shiftResidue = shiftResidue + shift;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    PPm_shift = PPm_shift + shift.';
    

    
    [beta, error_shift] = regressionFourthOrder(PPm_shift);
    betaError = [betaError norm(beta(4:15),1)];
    
    % 익명 함수 정의
    f2 = beta(1)*x.^4 + beta(2)*y.^4 + beta(3)*z.^4 + ...
         beta(4)*x.^2.*y.^2 + beta(5)*x.^2.*z.^2 + beta(6)*y.^2.*z.^2 + ...
         beta(7)*(x.^3).*y + beta(8)*(x.^3).*z + beta(9)*(y.^3).*x + ...
         beta(10)*(y.^3).*z + beta(11)*(z.^3).*x + beta(12)*(z.^3).*y - 1 + ...
         beta(13)*(x.^2).*y.*z + beta(14)*(y.^2).*z.*x + beta(15)*(z.^2).*x.*y;
    
    %소요시간 0.034초 내외
    %{
    figure;
    hold on;
    fimplicit3(f, [-1.5 1.5 -1.5 1.5 -1.5 1.5]);
    scatter3(PPm_shift(:,1), PPm_shift(:,2), PPm_shift(:,3), colors{mod(i,3)+1});
    hold off;
    axis equal;
    %}
end


%{
figure
hold on
f2_numeric = matlabFunction(f2,'Vars',[x,y,z]);
fimplicit3(f2_numeric,[-3.5 3.5 -3.5 3.5 -3.5 3.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),'b');
hold off
axis equal
%}
figure
plot(shiftResidueSet(1,:))
figure
plot(shiftResidueSet(2,:))