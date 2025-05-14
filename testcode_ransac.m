PP01 = generateRandomPointsOnHexagonPrism(7800)+randn(7800,3)*0.005;
PP02 = generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 = generateRandomPointsOnCylinder(8001)+randn(8001,3)*0.003;
PP04 = generateRandomPointsOnCylinder(2200)+randn(2200,3)*0.02;
PP05 = generateRandomPointsOnIcosahedron(4000);


t = linspace(0, 5400, 2000).';
Orbit = [ t(:), zeros(2000,3)];
Orbit(1,(2:4)) = [0.3 0.1 0];
w = 0.06;
n = [0 sin(35/57.92) cos(35/57.92)];
quater0= [cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];


%%
syms x y z
order = 6;
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



%% field of objective function
PPm_body = PP01;   
PPm_body(:,1) = PPm_body(:,1) * 1;
PPm_body(:,2) = PPm_body(:,2) * 1;
PPm_body(:,3) = PPm_body(:,3 )* 1.5;


rotm= quat2rotm([1,0,0,0]);
PPmR_body = PPm_body * rotm.';
PPm_use = PPmR_body - [Orbit(1,2), Orbit(1,3),0];

PPm_body2 = PP04;
PPm_body2(:,3) = PPm_body2(:,3 )* 4.5;
rotm= quat2rotm([1/sqrt(2),1/sqrt(2),0,0]);
PPmR_body2 = PPm_body2 * rotm.'- [Orbit(1,2), Orbit(1,3),0];
PPm_use = [PPm_use; PPmR_body2 ];


[beta_values,error]  = regressionFourthOrder(PPm_use,FuncsB);
inlierMask = abs(error) < 0.1;
score      = sum(inlierMask);
f1 = TermsB * beta_values;
[coeffsf0, monomialf0] = coeffs(f1 , [x,y,z]);   
Binit = double(coeffsf0);


ransacPar = struct('maxIter',600,'conf',0.90,'thresh',0.10,'minInlierRatio',0.6);

residuals1 = abs( FuncsB(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3))*beta_values-1); % N × #term
inlierMask1 = residuals1 < ransacPar.thresh;
score1      = sum(inlierMask1);





[DispRAN,BetaRAN, inlierMaskRAN] = PoliNavigationSolver3_Ransac(0,PPm_use,order,ransacPar);

PPm_use_unbias = PPm_use - DispRAN;
residuals2 = abs( FuncsB(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3))*BetaRAN-1); % N × #term
inlierMask2 = residuals2 < ransacPar.thresh;
score2      = sum(inlierMask2);

inlierMask3 = residuals2 < 0.1;
score3      = sum(inlierMask3);

%%
figure(1)
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),3,PPm_use(:,3),'filled');
hold on
fimplicit3(f1-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])



%%
f2_RAN = TermsB * BetaRAN;
PP_inlier = PPm_use(logical(inlierMaskRAN), :);
PP_inlier_unbias = PP_inlier - DispRAN;
figure(2)
scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,PP_inlier_unbias(:,3),'filled');
hold on
fimplicit3(f2_RAN-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

figure(3)
scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),3,PPm_use_unbias(:,3),'filled');
hold on
%fimplicit3(f2_RAN-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])
%%
colors = [0.2 0.6 0.2];   % FileExchange 함수
xr = linspace(-2.5, 2.5, 20);
yr = linspace(-2.5, 2.5, 20);
zr = linspace(-2.5, 2.5, 20);
[xg, yg, zg] = ndgrid(xr, yr, zr);
xi = xg(:) - DispRAN(1);
yi = yg(:) - DispRAN(2);
zi = zg(:) - DispRAN(3);

F   = FuncsB(xi,yi,zi) * BetaRAN - 1;
Vol = reshape(F, size(xg)); 

figure(4)
axis equal;
view([1 1 1]);
fv = isosurface(xg, yg, zg, Vol, 0);

patchH(m) = patch( fv, ...
        'FaceColor', colors(m,:), ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', 'none');
hold on 
scatter3(PPm_use(:,1), PPm_use(:,2), PPm_use(:,3), ...
         5, [0.8 0.8 0.8], '.');
hold off

%%
for m = 1:1
    % ① 좌표계 맞추기
    d   = DispRAN;               % 1×3
    XX  = X-d(1);  YY = Y-d(2);  ZZ = Z-d(3);

    % ② 다항식 값
    F   = FuncsB(XX,YY,ZZ) * BetaRAN - 1;

    % ③ 면 patch
    p   = patch( isosurface(X,Y,Z, abs(F), 0) );
    set(p, 'FaceColor', colors(m,:), ...
           'FaceAlpha', 0.15, ...
           'EdgeColor', 'none');
end
% ④ 데이터 점 (희미하게)
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),2,[0.5 0.5 0.5],'filled');
axis equal; view(3); camlight; lighting gouraud;
hold off

