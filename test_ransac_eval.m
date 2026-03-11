PP01 = generateRandomPointsOnHexagonPrism(7800)+randn(7800,3)*0.007;
PP02 = generateRandomPointsOnSpaceship(4000)+randn(4000,3)*0.006;
PP03 = generateRandomPointsOnCylinder(8001)+randn(8001,3)*0.005;
PP04 = generateRandomPointsOnCylinder(500)+randn(500,3)*0.02;
PP05 = generateRandomPointsOnIcosahedron(4000);

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P2 = P(P(:,2) >= 2 & P(:,2) <= 7 & P(:,3) >= 0 & P(:,3) <= 5,:);

t = linspace(0, 5400, 2000).';
Orbit = [ t(:), zeros(2000,3)];
Orbit(1,(2:4)) = [0.1 0.5 0];
w = 0.06;
n = [0 sin(35/57.92) cos(35/57.92)];
quater0= [cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];


%%
syms x y z 
order = 2;
%TermsC = homogeneFischerTerms(order);
%FuncsB = matlabFunction(TermsB);
TermsC = homogeneFischerTerms(order);
FuncsC = matlabFunction(TermsC);
Qinit = [1;0.0;0.0;0.00];
N = length(TermsC);


%% field of objective function
PPm_body = P2;   
PPm_body(:,1) = PPm_body(:,1) * 1;
PPm_body(:,2) = PPm_body(:,2) * 1;
PPm_body(:,3) = PPm_body(:,3 )* 2;


rotm= quat2rotm([1,0,0,0]);
PPmR_body = PPm_body * rotm.';
PPm_use = PPmR_body - [Orbit(1,2), Orbit(1,3),Orbit(1,4)];

PPm_body2 = PP04;
PPm_body2(:,3) = PPm_body2(:,3 )* 2.5;
rotm= quat2rotm([2/sqrt(9),2/sqrt(9),0,1/sqrt(9)]);
PPmR_body2 = PPm_body2 * rotm.' + [Orbit(1,3), Orbit(1,4), Orbit(1,2)];
%PPm_use = [PPm_use; PPmR_body2 ];


[beta_values,error]  = regressionFourthOrder(PPm_use,FuncsC);
inlierMask = abs(error) < 0.1;
score      = sum(inlierMask);
f1 = TermsC * beta_values;
[coeffsf0, monomialf0] = coeffs(f1 , [x,y,z]);   
Binit = double(coeffsf0);


%%

f1_RAN = TermsC * beta_values;
figure(1)
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),1,'b','filled');
hold on
fimplicit3(f1_RAN-1,[-100 100 -100 100 -100 100],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 10.5])



%% 구식 ransac 
%{

ransacPar = struct('maxIter',20600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.61);

residuals1 = abs( FuncsC(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3))*beta_values-1); % N × #term
inlierMask1 = residuals1 < ransacPar.thresh;
score1      = sum(inlierMask1);


center_shift = mean(PPm_use,1);
PPm_use_uncenter = PPm_use - center_shift;

[DispRAN,BetaRAN, inlierMaskRAN] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,ransacPar);
scoreRAN      = sum(inlierMaskRAN);

if ~isempty(BetaRAN)
    PPm_use_unbias = PPm_use - DispRAN;
    residuals2 = abs( FuncsC(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3))*BetaRAN-1); % N × #term
    inlierMask2 = residuals2 < ransacPar.thresh;
    score2      = sum(inlierMask2);
end


%}
%% 신형 ransac _ 속도 증가 by funcs 


ransacPar = struct('maxIter',20600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.61,'locIters',7, 'damping',0.85,'reg',1e-6);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);


residuals1 = abs( FuncsC(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3))*beta_values-1); % N × #term
inlierMask1 = residuals1 < ransacPar.thresh;
score1      = sum(inlierMask1);


center_shift = mean(PPm_use,1);
PPm_use_uncenter = PPm_use - center_shift;
%시행
[DispRAN,BetaRAN, inlierMaskRAN] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN      = sum(inlierMaskRAN);

if ~isempty(BetaRAN)
    PPm_use_unbias = PPm_use - DispRAN;
    residuals2 = abs( FuncsC(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3))*BetaRAN-1); % N × #term
    inlierMask2 = residuals2 < ransacPar.thresh;
    score2      = sum(inlierMask2);
end


%}
%%
if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN);
    PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = Funcs1_RAN(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,[-100 100 -100 100 -100 100],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-30.5 30.5 -10.5 10.5 -10.5 10.5])

    
    figure(3)
    scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,'b','filled');
    fimplicit3(f1_RAN-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor','none','FaceAlpha',0.5);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-30.5 30.5 -10.5 10.5 -10.5 10.5])

end
%%
colors = lines(10);   % FileExchange 함수
xr = linspace(-2.5, 2.5, 20);
yr = linspace(-2.5, 2.5, 20);
zr = linspace(-2.5, 2.5, 20);
[xg, yg, zg] = ndgrid(xr, yr, zr);
xi = xg(:) - DispRAN(1);
yi = yg(:) - DispRAN(2);
zi = zg(:) - DispRAN(3);

F   = FuncsB(xi,yi,zi) * BetaRAN - 1;
Vol = reshape(F, size(xg)); 

figure(6)
axis equal;
view([1 1 1]);
fv = isosurface(xg, yg, zg, Vol, 0);

patchH(m) = patch( fv, ...
        'FaceColor', colors(1,:), ...
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

