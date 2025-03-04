%2024 08 20
PP01 =generateRandomPointsOnHexagonPrism(2000)+randn(2000,3)*0.0001;
PP02 =generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 =generateRandomPointsOnCylinder(200)+randn(200,3)*0.0005;

shiftReal =  [0.0 -0 -0.01];
PPm_body = PP01 ;
PPm_body(:,3) = PPm_body(:,3) * 1.6;
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

q= [1 0 0 0];
%q= [cos(25/57.92) 0 sin(25/57.92) 0];
rotm= quat2rotm(q);
PPmR_body = PPm_body * rotm.';


figure
hold on
h = PPm_body(:,3);
scatter3(PPm_body(:,1),PPm_body(:,2),PPm_body(:,3),1,h,'filled');
h = PPm_panel1(:,3);
scatter3(PPm_panel1(:,1),PPm_panel1(:,2),PPm_panel1(:,3),1,h,'filled');
h = PPm_panel1(:,3);
scatter3(PPm_panel2(:,1),PPm_panel2(:,2),PPm_panel2(:,3),1,h,'filled');
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
view([1,1,1])
xlim([-6 6]);  
ylim([-6 6]);
zlim([-4 4]);
%title('Pointcloud of KOMPSAT-1 Model')


E1=load("equation_body.mat",'f3_translated_expanded');
E2=load("equation_panel1.mat","f3_translated_expanded");
E3=load("equation_panel2.mat","f3_translated_expanded");

figure
hold on
h = PPm_body(:,3);
scatter3(PPm_body(:,1),PPm_body(:,2),PPm_body(:,3),1,h,'filled');
h = PPm_panel1(:,3);
scatter3(PPm_panel1(:,1),PPm_panel1(:,2),PPm_panel1(:,3),1,h,'filled');
h = PPm_panel1(:,3);
scatter3(PPm_panel2(:,1),PPm_panel2(:,2),PPm_panel2(:,3),1,h,'filled');
fimplicit3(E1.f3_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
fimplicit3(E2.f3_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.05, 0.2, 0.5]);
fimplicit3(E3.f3_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.05, 0.2, 0.5]);
hold off


colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
view([1 1 1]);



%% 비동차항 먼저
syms x;
syms y;
syms z;
order = 6;
TermsA = nonhomogeneTerm(order);
betaA = sym('beta', [1, length(TermsA)]);
f1 = sum(betaA .* TermsA);

%점 입력
PPm_use = PPm_body;
[beta_values0,error0]    = regressionFourthOrder(PPm_use,TermsA);
f1_numeric = matlabFunction(f1, 'Vars', {[x, y, z], betaA});
f1_partial = @(x,y,z) f1_numeric([x,y,z], beta_values0.')-1;

%f1_total   = subs(f1,betaA,beta_values0);
%syms a b c 
%f1_t_shift = subs(f1_total, [x, y, z], [x-a, y-b, z-c])

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
%% 그냥 좌표 평균 찾기 겸 점 초기세팅

center_shift = median(PPm_use,1);
PPm_shift = PPm_use - center_shift;

shiftResidue = shiftReal.'
shiftResidue = shiftResidue - center_shift.'

%% Iteration 초기 설정
%shifted_f = subs(poly, [x, y, z], [x - a, y - b, z - c])
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

[beta_values,error]    = regressionFourthOrder(PPm_shift,TermsB);

error_shift = error;  % 초기 에러 설정

shiftSet = [];
shift1Set = [];
shiftResidueSet = [];
shiftSum = [0 ; 0 ; 0];
numIterations = 10;  % 반복 횟수 설정
shiftSum = center_shift.';


%% 그림 그리기_initial  
f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
f2_substituted = subs(f2, betaB, beta_values.');
f2_translated = subs(f2_substituted-1, [x,y,z], [x-center_shift(1), y-center_shift(2), z-center_shift(3)]);
f2_translated_expanded = expand(f2_translated);
f2_handle = matlabFunction(f2_translated_expanded, 'Vars', [x,y,z]);
[coeffs1, monomial] = coeffs(f2_partial , [x,y,z])
%coeffs(f2_substituted , [x,y,z])

figure
h2 = scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),1,PPm_shift(:,3),'filled');
colormap(jet);
hold on
h1 = fimplicit3(f2_partial,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
hold off

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
%title('Pointcloud Model - Panel')
%title('Pointcloud Model - Body')

figure
h2 = scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),1,PPm_use(:,3),'filled');
colormap(jet);
hold on
h1 = fimplicit3(f2_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.05, 0.2, 0.5]);
hold off

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
%title('Pointcloud Model - Panel')
%{
figure
h2 = scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),1,PPm_shift(:,3),'filled');
colormap(jet);
hold on
f2_ngd = @(x,y,z) 1.10736*x^6 + 2.46766*y^6+ 0.0163682*x^5*y + 12.3976*x^4*y^2 - 0.434712*x^3*y^3 - 2.58541*x^2*y^4  + 0.143499*x*y^5 - 0.286339*x^4*y*z  - 0.0736799*x*y^4*z - 0.0199605*y^5*z - 0.11549*x^4*z^2 - 0.254459*x^2*y^2*z^2  - 0.210524*y^4*z^2 - 0.093437*x^2*z^4 + 0.016984*x*y*z^4 - 0.0840934*y^2*z^4  + 0.0635403*z^6-1;
h1 = fimplicit3(f2_ngd,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
hold off
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
%}
%%
f2_2D=  @(x,y) -f2_partial(x,y,0);
[X, Y] = meshgrid(linspace(-2.5, 2.5, 100), ...
                  linspace(-2.5, 2.5, 100));
figure
fcontour(f2_2D,[-1.5 1.5 -1.5 1.5],'fill','on','LevelStep',0.2);
%caxis auto;
caxis ([-1 1]);
colormap(jet(256));
xlabel('$\mathrm{X}$', 'Interpreter', 'latex')
ylabel('$\mathrm{Y}$', 'Interpreter', 'latex')
figure
fsurf(f2_2D,[-1.5 1.5 -1.5 1.5])
caxis([-1.5 1]);
zlim([ -1.2 1.1]);
colormap("jet");
xlabel('$\mathrm{X}$', 'Interpreter', 'latex')
ylabel('$\mathrm{Y}$', 'Interpreter', 'latex')
zlabel('$\mathrm{g(X,Y)}$', 'Interpreter', 'latex')
view([1,1,1])

f2_2Dcut = @(x,y) max(f2_2D(x,y),0);

figure
fsurf(f2_2Dcut,[-1.5 1.5 -1.5 1.5])
caxis([-0.5 1]);
zlim([ -0.2 1.1]);
xlabel('$\mathrm{X}$', 'Interpreter', 'latex')
ylabel('$\mathrm{Y}$', 'Interpreter', 'latex')
zlabel('$\mathrm{ReLU(g(X,Y))}$', 'Interpreter', 'latex')
colormap("jet");

f2_2Dsoftplus = @(x,y) log(exp(f2_2D(x,y)*3)+1)/3;
figure
fsurf(f2_2Dsoftplus,[-1.5 1.5 -1.5 1.5])
caxis([-0.5 1]);
zlim([ -0.2 1.1]);
xlabel('$\mathrm{X}$', 'Interpreter', 'latex')
ylabel('$\mathrm{Y}$', 'Interpreter', 'latex')
zlabel('$\mathrm{ReLU(g(X,Y))}$', 'Interpreter', 'latex')
colormap("jet");

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
    shiftSum = shiftSum + shift ;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    
    PPm_shift = PPm_shift + shift.';
    [beta_values, error_shift] = regressionFourthOrder(PPm_shift,TermsB);

    %0.001s 
    
end

f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
figure
hold on
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3));
%scatter3(PP02(:,1),PP02(:,2),PP02(:,3),'b');
hold off
colormap(jet);
axis equal
xlim([-4 4]);
ylim([-4 4]);
title('After')

f3_substituted = subs(f2, betaB, beta_values.');
f3_translated = subs(f3_substituted-1, [x,y,z], [x-shiftSum(1), y-shiftSum(2), z-shiftSum(3)]);
f3_translated_expanded = expand(f3_translated);


figure
h2 = scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),1,PPm_use(:,3),'filled');
colormap(jet);
hold on
h1 = fimplicit3(f3_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.05, 0.2, 0.5]);
hold off

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal

%{
norm(error_shift,1)
evaluateModel(PPm_shift,f2_numeric, beta_values)
%}
figure
subplot(1,3,1)
plot(shiftResidueSet(1,:))
subplot(1,3,2)
plot(shiftResidueSet(2,:))
subplot(1,3,3)
plot(shiftResidueSet(3,:))
