PP01 = generateRandomPointsOnHexagonPrism(7800)+randn(7800,3)*0.007;
PP02 = generateRandomPointsOnSpaceship(4000)+randn(4000,3)*0.006;
PP03 = generateRandomPointsOnCylinder(8001)+randn(8001,3)*0.005;
PP04 = generateRandomPointsOnCylinder(500)+randn(500,3)*0.02;
PP05 = generateRandomPointsOnIcosahedron(4000);

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8,:);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 &  P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13,:);



figure(11)
scatter3(P1(:,1),P1(:,2),P1(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 15.5 -25.5 20.5])



t = linspace(0, 5400, 2000).';
Orbit = [ t(:), zeros(2000,3)];
Orbit(1,(2:4)) = [0.1 0.5 0];
w = 0.06;
n = [0 sin(35/57.92) cos(35/57.92)];
quater0= [cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];


%%
syms x y z 
order = 4;
%TermsC = homogeneFischerTerms(order);
%FuncsB = matlabFunction(TermsB);
TermsC = homogeneFischerTerms(order);
%TermsC2 = homogeneFischerTerms(order);
FuncsC = matlabFunction(TermsC);
Qinit = [1;0.0;0.0;0.00];
N = length(TermsC);
%N2 = length(TermsC);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);
warning('off','MATLAB:nearlySingularMatrix')
%% field of objective function
PPm_body = P4;   
PPm_body(:,1) = PPm_body(:,1) * 1;
PPm_body(:,2) = PPm_body(:,2) * 1;
PPm_body(:,3) = PPm_body(:,3 )* 1;


rotm= quat2rotm([1,0,0,0]);
PPmR_body = PPm_body * rotm.';
PPm_use = PPmR_body - [Orbit(1,2), Orbit(1,3),Orbit(1,4)];

[beta_values,error]  = regressionFourthOrder(PPm_use,FuncsC);
inlierMask = abs(error) < 0.1;
score      = sum(inlierMask);
f1 = TermsC * beta_values;
[coeffsf0, monomialf0] = coeffs(f1 , [x,y,z]);   
Binit = double(coeffsf0);


%%
%{
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
axis([-30.5 30.5 -10.5 10.5 -10.5 15.5])
%}
%% 신형 ransac _ 속도 증가 by funcs 


ransacPar = struct('maxIter',7600,'conf',0.96,'thresh',0.554,'minInlierRatio',0.65,'updateThresh',0.62,'locIters',7, 'damping',0.85,'reg',1e-6);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);


residuals1 = abs( FuncsC(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3))*beta_values-1); % N × #term
inlierMask1 = residuals1 < ransacPar.thresh;
score1      = sum(inlierMask1);


%center_shift = mean(PPm_use,1);
PPm_use_uncenter = PPm_use;
%시행
warning('off','MATLAB:nearlySingularMatrix')
[DispRAN4,BetaRAN4, inlierMaskRAN4] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN4      = sum(inlierMaskRAN4);

if ~isempty(BetaRAN4)
    PPm_use_unbias = PPm_use - DispRAN4;
    residuals2 = abs( FuncsC(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3))*BetaRAN4-1); % N × #term
    inlierMask2 = residuals2 < ransacPar.thresh;
    score4      = sum(inlierMask2);
end


%%

figure(11)
scatter3(P1(:,1),P1(:,2),P1(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-60.5 60.5 -35.5 45.5 -20.5 30.5])

PPm_use = P1;  
ransacPar = struct('maxIter',9600,'conf',0.96,'thresh',0.714,'minInlierRatio',0.65,'updateThresh',0.69,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN1,BetaRAN1, inlierMaskRAN1] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN1      = sum(inlierMaskRAN1);

if ~isempty(BetaRAN1)
    PPm_use_unbias = PPm_use - DispRAN1;

end
%%

figure(11)
scatter3(P2(:,1),P2(:,2),P1(:,3),3,'r','filled');
hold on 
scatter3(P3(:,1),P3(:,2),P3(:,3),3,'r','filled');
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P2;  
ransacPar = struct('maxIter',8600,'conf',0.96,'thresh',0.734,'minInlierRatio',0.65,'updateThresh',0.69,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN2,BetaRAN2, inlierMaskRAN2] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN2      = sum(inlierMaskRAN2);

if ~isempty(BetaRAN2)
    PPm_use_unbias = PPm_use - DispRAN2;

end
PPm_use = P3;  
ransacPar = struct('maxIter',20600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);
[DispRAN3,BetaRAN3, inlierMaskRAN3] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
scoreRAN3      = sum(inlierMaskRAN3);
if ~isempty(BetaRAN3)
    PPm_use_unbias = PPm_use - DispRAN3;
end
%%
P5 = P(P(:,1) >= -3.5 & P(:,1) <= 3.5 &  P(:,2) <= 2.5 & P(:,3) >= 11.5 & P(:,3) <= 18.5,:);
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.55,'locIters',7, 'damping',0.85,'reg',1e-6);

figure(11)
scatter3(P5(:,1),P5(:,2),P5(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P5;  

%시행
[DispRAN5,BetaRAN5, inlierMaskRAN5] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN5      = sum(inlierMaskRAN5);

if ~isempty(BetaRAN5)
    PPm_use_unbias = PPm_use - DispRAN5;
    residualsaf = abs( FuncsC(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3))*BetaRAN5-1); % N × #term
    inlierMaskaf = residualsaf < ransacPar.thresh;
    score5      = sum(inlierMaskaf);
end


%%

P6 = P(P(:,1) >= 2 & P(:,1) <= 13 &  P(:,2) <= 2.5 & P(:,3) >= 13 ,:);
%ransacPar = struct('maxIter',10600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.55,'locIters',7, 'damping',0.85,'reg',1e-6);

figure(11)
scatter3(P6(:,1),P6(:,2),P6(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P6;  

%시행
[DispRAN6,BetaRAN6, inlierMaskRAN6] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
scoreRAN6      = sum(inlierMaskRAN6);

if ~isempty(BetaRAN6)
    PPm_use_unbias = PPm_use - DispRAN6;

end


%%
P7 = P(P(:,1) >= -10 & P(:,1) <=-2 &  P(:,2) <= 2.5 & P(:,3) >= 13 ,:);

figure(11)
scatter3(P7(:,1),P7(:,2),P7(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P7;  
[DispRAN7,BetaRAN7, inlierMaskRAN7] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN7      = sum(inlierMaskRAN7);

if ~isempty(BetaRAN7)
    PPm_use_unbias = PPm_use - DispRAN7;

end
%%
P8 =  P(P(:,1) >= -3.0 & P(:,1) <= 3.0 &  P(:,2) <= 2.6 & P(:,2) >= -2.5 & P(:,3) >= -3 & P(:,3) <= 3,:);
ransacPar = struct('maxIter',13600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.55,'locIters',7, 'damping',0.85,'reg',1e-6);

figure(11)
scatter3(P8(:,1),P8(:,2),P8(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P8;  
[DispRAN8,BetaRAN8, inlierMaskRAN8] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN8      = sum(inlierMaskRAN8);

if ~isempty(BetaRAN8)
    PPm_use_unbias = PPm_use - DispRAN8;

end
%%
P9 =  P(P(:,1) >= -6.5 & P(:,1) <= -2.5 &  P(:,2) <= 3 & P(:,2) >= -3 & P(:,3) >= -3 & P(:,3) <= 3,:);

figure(11)
scatter3(P9(:,1),P9(:,2),P9(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P9;  
ransacPar = struct('maxIter',20600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.55,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN9,BetaRAN9, inlierMaskRAN9] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN9      = sum(inlierMaskRAN9);

if ~isempty(BetaRAN9)
    PPm_use_unbias = PPm_use - DispRAN9;

end

%%
P10 =  P(P(:,1) >= 2.5 & P(:,1) <= 9 &  P(:,2) <= 3 & P(:,2) >= -3 & P(:,3) >= -3 & P(:,3) <= 3,:);

figure(11)
scatter3(P10(:,1),P10(:,2),P10(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P10;  
%ransacPar = struct('maxIter',12600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.55,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN10,BetaRAN10, inlierMaskRAN10] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN10      = sum(inlierMaskRAN10);

if ~isempty(BetaRAN10)
    PPm_use_unbias = PPm_use - DispRAN10;

end

%%
P11 =  P(P(:,1) >= -2.5 & P(:,1) <= 3 &  P(:,2) <= 7.3 & P(:,2) >= 2.5 & P(:,3) >= -3 & P(:,3) <= 3,:);

figure(11)
scatter3(P11(:,1),P11(:,2),P11(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P11;  
%ransacPar = struct('maxIter',12600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN11,BetaRAN11, inlierMaskRAN11] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN11      = sum(inlierMaskRAN11);

if ~isempty(BetaRAN11)
    PPm_use_unbias = PPm_use - DispRAN11;

end
%%
P12 =  P(P(:,1) >= -2.5 & P(:,1) <= 3 &  P(:,2) <= -2.5 & P(:,2) >= -9.5 & P(:,3) >= -3 & P(:,3) <= 3,:);

figure(11)
scatter3(P12(:,1),P12(:,2),P12(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -10.5 20.5])

PPm_use = P12;  
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN12,BetaRAN12, inlierMaskRAN12] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN12      = sum(inlierMaskRAN12);

if ~isempty(BetaRAN12)
    PPm_use_unbias = PPm_use - DispRAN12;

end

%%
P13 =  P(P(:,1) >= -2.5 & P(:,1) <= 3 &  P(:,2) <= 3.0 & P(:,2) >= -2.0 &  P(:,3) <= -5,:);
figure(11)
scatter3(P13(:,1),P13(:,2),P13(:,3),3,'r','filled');
hold on 
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-30.5 30.5 -10.5 10.5 -30.5 20.5])

PPm_use = P13;  
%ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN13,BetaRAN13, inlierMaskRAN13] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN13      = sum(inlierMaskRAN13);

if ~isempty(BetaRAN13)
    PPm_use_unbias = PPm_use - DispRAN13;

end

%% end wing -> not feasible , manual oper required
P14 =  P( P(:,1) >= 27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)>=13,:);
P15 =  P( P(:,1) >= 27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)<=-7,:);
P16 =  P( P(:,1) <= -27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)>=13,:);
P17 =  P( P(:,1) <= -27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)<=-7,:);
figure(11)
scatter3(P14(:,1),P14(:,2),P14(:,3),3,'r','filled');
hold on 
scatter3(P15(:,1),P15(:,2),P15(:,3),3,'b','filled');
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-40.5 40.5 -30.5 30.5 -30.5 20.5])

PPm_use = P16;  
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.554,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN16,BetaRAN16, inlierMaskRAN16] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN16     = sum(inlierMaskRAN16);

if ~isempty(BetaRAN16)
    PPm_use_unbias = PPm_use - DispRAN16;

end

%% side bottom panel
P18 =  P( P(:,1) >= 27 & P(:,1) <= 35 & P(:,2)+1.8*P(:,3)<=8.5 ,:);
P19 =  P( P(:,1) >= 40 & P(:,1) <= 50 & P(:,2)+1.8*P(:,3)<=8.5 ,:);

P20 =  P( P(:,1) <= -27 & P(:,1) >= -35 & P(:,2)+1.8*P(:,3)<=8.5 ,:);
P21 =  P( P(:,1) <= -40 & P(:,1) >= -50 & P(:,2)+1.8*P(:,3)<=8.5 ,:);
figure(11)
scatter3(P18(:,1),P18(:,2),P18(:,3),3,'r','filled');
hold on 
scatter3(P19(:,1),P19(:,2),P19(:,3),3,'b','filled');
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-50.5 50.5 -30.5 30.5 -30.5 20.5])

PPm_use = P18;  
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.754,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN18,BetaRAN18, inlierMaskRAN18] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN18      = sum(inlierMaskRAN18);

if ~isempty(BetaRAN18)
    PPm_use_unbias = PPm_use - DispRAN18;

end
%%
PPm_use = P19;  
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.704,'minInlierRatio',0.65,'updateThresh',0.65,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN19,BetaRAN19, inlierMaskRAN19] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN19      = sum(inlierMaskRAN19);

if ~isempty(BetaRAN19)
    PPm_use_unbias = PPm_use - DispRAN19;

end
%%
PPm_use = P20;  
[DispRAN20,BetaRAN20, inlierMaskRAN20] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN20      = sum(inlierMaskRAN20);

if ~isempty(BetaRAN20)
    PPm_use_unbias = PPm_use - DispRAN20;

end
%%
PPm_use = P21;  
[DispRAN21,BetaRAN21, inlierMaskRAN21] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN21      = sum(inlierMaskRAN21);

if ~isempty(BetaRAN21)
    PPm_use_unbias = PPm_use - DispRAN21;

end
%% side bar
P22 =  P( P(:,1) >= 27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)<=8& 1.8*P(:,2)-P(:,3)>=-2,:);

P23 =  P( P(:,1) <= -27 & P(:,2)+2*P(:,3)>=6 & 1.8*P(:,2)-P(:,3)<=8& 1.8*P(:,2)-P(:,3)>=-3,:);

figure(11)
scatter3(P22(:,1),P22(:,2),P22(:,3),3,'r','filled');
hold on 
scatter3(P23(:,1),P23(:,2),P23(:,3),3,'b','filled');
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-40.5 40.5 -30.5 30.5 -30.5 20.5])

PPm_use = P22;  
ransacPar = struct('maxIter',8600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.75,'updateThresh',0.63,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN22,BetaRAN22, inlierMaskRAN22] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN22      = sum(inlierMaskRAN22);

if ~isempty(BetaRAN22)
    PPm_use_unbias = PPm_use - DispRAN22;

end
PPm_use = P23;  
ransacPar = struct('maxIter',8600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.75,'updateThresh',0.62,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN23,BetaRAN23, inlierMaskRAN23] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN23      = sum(inlierMaskRAN23);

if ~isempty(BetaRAN23)
    PPm_use_unbias = PPm_use - DispRAN23;

end
%%
P24 =  P( P(:,1) <= 27 &P(:,1) >= 23 & P(:,3) <10 ,:);
P25 =  P( P(:,1) >= -27 & P(:,1) <= -20 &P(:,3) <8,:);
figure(11)
scatter3(P24(:,1),P24(:,2),P24(:,3),3,'r','filled');
hold on 
scatter3(P25(:,1),P25(:,2),P25(:,3),3,'b','filled');
scatter3(P(:,1),P(:,2),P(:,3),1,'black','filled');
hold off
axis equal
axis([-40.5 40.5 -30.5 30.5 -30.5 20.5])

PPm_use = P24;  
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.80,'updateThresh',0.68,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN24,BetaRAN24, inlierMaskRAN24] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN24      = sum(inlierMaskRAN24);

if ~isempty(BetaRAN24)
    PPm_use_unbias = PPm_use - DispRAN24;

end

PPm_use = P25;  
%ransacPar = struct('maxIter',8600,'conf',0.96,'thresh',0.454,'minInlierRatio',0.55,'updateThresh',0.68,'locIters',7, 'damping',0.85,'reg',1e-6);

[DispRAN25,BetaRAN25, inlierMaskRAN25] = PoliNavigationSolver3_FischerRansac(0,PPm_use,order,length(TermsC),Funcs1,Grads1,ransacPar);
%
scoreRAN25      = sum(inlierMaskRAN25);

if ~isempty(BetaRAN25)
    PPm_use_unbias = PPm_use - DispRAN25;

end

%% graph for any ransac

%target setting (whether if 4 or 5)
BetaRAN = BetaRAN8;
DispRAN = DispRAN8;
inlierMaskRAN = inlierMaskRAN8;


if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN);

    PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = Funcs1_RAN(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,[-20 20 -10 10 -30 30],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 150);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-25.5 25.5 -10.5 10.5 -15.5 15.5])

    
    figure(3)
    scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,'b','filled');
    
    fimplicit3(f1_RAN-1,[-20 20 -20 20 -30 30],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor','none','FaceAlpha',0.5, 'MeshDensity', 150);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-25.5 25.5 -15.5 15.5 -15.5 15.5])

end


%% all plot
figure(12)
clf
% 배경 점 (전체 점군)
scatter3(P(:,1), P(:,2), P(:,3), 1, 'k', 'filled');
hold on

% 색상 팔레트 정의 (필요 시 반복 사용)
colors = lines(140);  % MATLAB 기본 컬러맵 중 하나

% 점군 리스트 (P2~P14 존재 가정)
for i = 4:8
    Pi = eval(sprintf('P%d', i));   % P2, P3, ..., P14 변수 참조
    scatter3(Pi(:,1), Pi(:,2), Pi(:,3), 3, colors(round(i*1.7),:), 'filled');
end

hold off
axis equal
axis([-40.5 40.5 -20.5 25.5 -30.5 30.5])
xlabel('X'); ylabel('Y'); zlabel('Z');
title('ISS Stationary Points with P2–P14 Overlays');

%% data save
%사실 여기에는 별거 없음. 이렇게 저장하면 큰일남. 
% ISS plot manual fix 갈것
%{
Data3(25) = struct('P',[],'Disp',[],'Beta',[],'InlierMask',[]);

for i = 1:25
    vB = sprintf('BetaRAN%d', i);
    if ~exist(vB,'var')
        fprintf('skip %d (Beta missing)\n', i); continue;
    end
    Beta = eval(vB);
    if isempty(Beta)
        fprintf('skip %d (Beta empty)\n', i); continue;
    end

    vD = sprintf('DispRAN%d', i);
    if ~exist(vD,'var'), fprintf('skip %d (Disp missing)\n', i); continue; end
    Disp = eval(vD);  Disp = Disp(:);  % [dx dy dz] 행벡터 보장
    vP = sprintf('P%d', i);
    if exist(vP,'var'), P = eval(vP); else, P = []; end
    vM = sprintf('inlierMaskRAN%d', i);
    if exist(vM,'var'), M = eval(vM); else, M = []; end

    Data3(i).P          = P;
    Data3(i).Disp       = Disp;
    Data3(i).Beta       = Beta;
    Data3(i).InlierMask = M;
    f_raw_sym = TermsC * Beta;
    f_raw     = matlabFunction(f_raw_sym, 'Vars', [x y z]);
    f_shifted = @(X,Y,Z) f_raw(X - Disp(1), Y - Disp(2), Z - Disp(3));
    
    fname = sprintf('%s/Funcs_ISS%d.py', outDir, i);
    fprintf('Exporting %s ...\n', fname);
    exportFuncToPython({f_shifted},fname, 1, {}, {'x','y','z'});
   
end

%}
%% panel 후처리


PPm_use = P17;
BetaRAN = BetaRAN2;
DispRAN = [-42.6,-14.2,16.7];
DispRAN17 = [-42.6,-14.2,16.7];
%[42.6,23.71,-4.7];%DispRAN14;
inlierMaskRAN = inlierMaskRAN17;


theta = pi/3;       % 60도
c = cos(theta);     % 0.5
s = sin(theta);     % 
R  = [1  0   0;0  c  -s;0  s   c]*[0  -1  0;1  0  0;0  0  1];

s_x = 4;   % 길이배수 (2면 x축 길이 2배)
s_y = 2.5;   % y축 스케일 (변경없음)
s_z = 1.8;   % z축 스케일 (변경없음)
S_inv = diag([1/s_x, 1/s_y, 1/s_z]);

Rt = R.';
M = S_inv * Rt;

if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN);

    %f1_shifted = @(X,Y,Z) Funcs1_RAN( X/s_x, Y/s_y, Z/s_z);
    f1_rotated = @(X,Y,Z) Funcs1_RAN( ...
    M(1,1).*X + M(1,2).*Y + M(1,3).*Z, ...
    M(2,1).*X + M(2,2).*Y + M(2,3).*Z, ...
    M(3,1).*X + M(3,2).*Y + M(3,3).*Z );


    PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = f1_rotated(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(@(x,y,z) f1_rotated(x,y,z) - 1,[-20 20 -20 20 -16 16],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 250);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 1.0]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-40.5 40.5 -20.5 25.5 -30.5 30.5])


end
%%
Data4(:) = Data3(:);
for i = 14:17

    vD = sprintf('DispRAN%d', i);
    if ~exist(vD,'var'), fprintf('skip %d (Disp missing)\n', i); continue; end
    Disp = eval(vD);  Disp = Disp(:);  % [dx dy dz] 행벡터 보장
    vP = sprintf('P%d', i);
    if exist(vP,'var'), P = eval(vP); else, P = []; end


    Data3(i).P          = P;
    Data3(i).Disp       = Disp;
    
    orig_x = M(1,1)*(x - Disp(1)) + M(1,2)*(y - Disp(2)) + M(1,3)*(z - Disp(3));
    orig_y = M(2,1)*(x - Disp(1)) + M(2,2)*(y - Disp(2)) + M(2,3)*(z - Disp(3));
    orig_z = M(3,1)*(x - Disp(1)) + M(3,2)*(y - Disp(2)) + M(3,3)*(z - Disp(3));

    f_sym = TermsC * BetaRAN2;
    f_rotated_sym = subs(f_sym, [x,y,z], [orig_x, orig_y, orig_z]);
    f_rotated_raw     = matlabFunction(f_rotated_sym, 'Vars', [x y z]);
    %f_ = @(X,Y,Z) f_raw( ...
    %M(1,1).*X + M(1,2).*Y + M(1,3).*Z, ...
    %M(2,1).*X + M(2,2).*Y + M(2,3).*Z, ...
    %M(3,1).*X + M(3,2).*Y + M(3,3).*Z );
    
    fname = sprintf('image_RANSAC/Funcs_ISS%d.py', i);
    fprintf('Exporting %s ...\n', fname);
    exportFuncToPython(f_rotated_raw,fname, {'x','y','z'});
   
end