

% (2) 시간 범위 설정 
t = linspace(0, 5400, 2000);
Orbit = [ t(:), zeros(2000,3)];


%%
PP01 = generateRandomPointsOnHexagonPrism(1000)+randn(1000,3)*0.02;
PP02 = generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 = generateRandomPointsOnCylinder(7000)+randn(7000,3)*0.0005;
PP04 = generateRandomPointsOnCube(200);

PPm_body = PP01 ;
PPm_body(:,1) = PPm_body(:,1) * 1.276;
PPm_body(:,2) = PPm_body(:,2) * 1.14;
PPm_body(:,3) = PPm_body(:,3) * 1.65;

order = 6;
quater0 = zeros(8000, 4);
quater0(:,1) = 1;
quaterinit=zeros(8000, 4);
for i = 1:8000
    qRand = randn(1, 4); 
    qRand = qRand / norm(qRand);
    
    quaterinit(i, :) = qRand;
end
quaterinit(2,:) = [1, 0, 0, 0]; 
%{
figure;
histogram(quater0(:, 2)); 
xlabel('q0');
ylabel('Count');
title('Histogram of q0');
%}
%%
syms x y z
TermsB = homogeneTerm(order);
FuncsB = matlabFunction(TermsB);
Qinit = [1;0.0;0.0;0.00];
N = length(TermsB);

syms q0 q1 q2 q3  real
q = [q0; q1; q2; q3];
% 회전행렬 R(q) 정의
Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
       2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
       2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];

% 치환된 변수 정의
xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;


betaR = sym('beta', [1, N]);
f0 = betaR * TermsB(:);
f0_rotated = subs(f0, [x, y, z], [xr, yr, zr]);

M_sym = sym(zeros(N));
TermsB_rot = subs(TermsB, [x,y,z],[xr,yr,zr]);
for j = 1:N
    [cList,mList] = coeffs(TermsB_rot(j),[x,y,z]);
    for k = 1:length(mList)
        idx = find(mList(k) == TermsB); 
        M_sym(idx,j) = cList(k);
    end
end
global Mhandle;
Mhandle = matlabFunction(M_sym, 'Vars', {q0,q1,q2,q3});



%%

[beta_values,error]  = regressionFourthOrder(PPm_body,FuncsB);
f1 = TermsB * beta_values;
[coeffsf0, monomialf0] = coeffs(f1 , [x,y,z]);   
Binit = double(coeffsf0);

w = ones(N,1);
w(1) = 100;
w(2) = 100;
w(3) = 100;
fObj1 = @(qVec) sum( w .* ((Mhandle(qVec(1),qVec(2),qVec(3),qVec(4))*beta_values - beta_values ).^2));
resultObj1 = arrayfun(@(idx) fObj1(quaterinit(idx,:)), 1:8000)';
figure(9)
scatter3(quaterinit(:,2),quaterinit(:,3),quaterinit(:,4),3,resultObj1(:),'filled');
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis equal
colorbar 
caxis ([0 300])

%%
range = 1:120;
initialFvals = zeros(120, 1);  % 초기 목적함수 값
finalFvals = zeros(120, 1);

Array_Disp=zeros(2000,4);
Array_QB=zeros(2000,5);
Array_BetaB=zeros(2000,N+1);

 
tic
for i = range
    %disp(i);
    PP01 = generateRandomPointsOnHexagonPrism(100)+randn(100,3)*0.02;
    PP02 = generateRandomPointsOnCube(100)+randn(100,3)*0.02;
    PPm_body = PP01 ;
    PPm_body(:,1) = PPm_body(:,1) * 1.276;
    PPm_body(:,2) = PPm_body(:,2) * 1.14;
    PPm_body(:,3) = PPm_body(:,3) * 1.65;


    PPmR_body = PPm_body;
    PP_use = PPmR_body - [Orbit(i,2), Orbit(i,3),0];
    %PP_use = PPmR_body;
    Qinit = quaterinit(i,:).';
    %disp(beta_values);
    [DispB,QB,BetaB,inValB,optValB] = PoliNavigationSolver2(0,PP_use,order,beta_values,Qinit);
    %disp(QB.');
    initialFvals(i) = inValB;
    finalFvals(i) = optValB;
    Array_Disp(i,:) = [Orbit(i,1), DispB];
    Array_QB(i,:) = [Orbit(i,1), QB.'];
    Array_BetaB(i,:) =[Orbit(i,1), BetaB];

end
toc
%% result
figure(1)
scatter3(PPm_body(:,1),PPm_body(:,2),PPm_body(:,3),3,PPm_body(:,3),'filled');
hold on
fimplicit3(f1-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
%%
%quaternion

figure(10)
scatter3(Array_QB(range,3),Array_QB(range,4),Array_QB(range,5),3,finalFvals(range),'filled');
hold on
scatter3(quater0(range,2),quater0(range,3),quater0(range,4),3,'k','filled');
hold off
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis equal
legend('Estimation','Ground Truth');
colorbar 
caxis ([0 100])
figure(11)
histogram(finalFvals(range), 'BinWidth', 10);
xlabel('finalFvals');
ylabel('Frequency');
title('Histogram of finalFvals');

qTemp = zeros(8000,4);
for i= range
    qTemp(i,:) = quatmultiply(quatconj(Array_QB(i,2:5)),quater0(i,:));
end
figure(12)
scatter3(qTemp(:,2),qTemp(:,3),qTemp(:,4),4,qTemp(:,1),'filled')
axis equal
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis ([-1 1 -1 1 -1 1]);


%%
TermsB = homogeneTerm(order);

for i = 84
    fB = Array_BetaB(i,2:N+1) *  TermsB(:);
    fB_shift =  subs(fB, [x, y, z], [x - Array_Disp(i,2), y - Array_Disp(i,3), z - Array_Disp(i,4)])-1;
    PP01 = generateRandomPointsOnHexagonPrism(100)+randn(100,3)*0.02;
    PP02 = generateRandomPointsOnCube(100)+randn(100,3)*0.02;
    PPm_body = PP01 ;
    PPm_body(:,1) = PPm_body(:,1) * 1.276;
    PPm_body(:,2) = PPm_body(:,2) * 1.14;
    PPm_body(:,3) = PPm_body(:,3) * 1.65;

    rotm= quat2rotm(quater0(i,:));
    PPmR_body = PPm_body * rotm.';
    PP_use = PPmR_body - [Orbit(i,2), Orbit(i,3),0];
end 

figure(4)
scatter3(PP_use(:,1),PP_use(:,2),PP_use(:,3),1,PP_use(:,3),'filled');
hold on
fimplicit3(fB_shift,[-30 30 -50 50 -20 30],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
grid on

%% rotated image
syms q0 q1 q2 q3 real
syms x y z
Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
       2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
       2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];

xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;

%fA_rot = subs(fA, [x,y,z], [xr,yr,zr]);
%fA_rot = subs(fA_rot, [q0,q1,q2,q3], QA.');
fB_rot = subs(fB, [x,y,z], [xr,yr,zr]);
fB_rot2 = subs(fB_rot, [q0,q1,q2,q3], Array_QB(i,2:5));
[coeffsfB, monomialfB] = coeffs(fB_rot2 , [x,y,z]);
coeffsfB = double(coeffsfB);
coeffsfB_unused = sym([]);


%RqA = double(subs(Rq,[q0,q1,q2,q3],[QA(1),-QA(2),-QA(3),-QA(4)]));
%PPmR_shiftA = PPmR_body - DispA;
%PPmR_rotA = (RqA * PPmR_shiftA.').';

RqB = double(subs(Rq,[q0,q1,q2,q3],[ Array_QB(i,2),-Array_QB(i,3),-Array_QB(i,4),- Array_QB(i,5)]));
PPmR_shiftB = PPmR_body;
PPmR_rotB = (RqB * PPmR_shiftB.').';


figure(3)
scatter3(PPmR_body(:,1),PPmR_body(:,2),PPmR_body(:,3),3,PPmR_body(:,3),'filled');
hold on
fimplicit3(fB-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])

figure(2)
scatter3(PPmR_rotB(:,1),PPmR_rotB(:,2),PPmR_rotB(:,3),3,PPmR_rotB(:,3),'filled');
hold on
fimplicit3(fB_rot2-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid on



