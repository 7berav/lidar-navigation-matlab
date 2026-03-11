%2024 08 20
%발표용 이미지 만들기 (대상 kompsat 몸통)

tic
N = 5012;

PP01 =generateRandomPointsOnCube(1403)+randn(1403,3)*0.0001;
PP12 =generateRandomPointsOnHexagonPrism(N);
PP13 =generateRandomPointsOnCylinder(N);
toc

PPim= PP01;
PPim(:,1)=PPim(:,1)*0.50;
PPim(:,2)=PPim(:,2)*0.5;
PPim(:,3)=PPim(:,3)*1;

shiftReal =  [0.0 -0.0 +0.00];
PPim = PPim + shiftReal;
%PPim(:,1)=PPim(:,1)*0.5;
%PPim(:,2)=PPim(:,2)*0.5;
%PPim(:,3)=PPim(:,3)*1.1;
%PPim = PPim+randn(N,3)*0.003;

%{
figure;

scatter3(PPim(:,1),PPim(:,2),PPim(:,3),'black');
axis equal
zlim([-1.5 1.5]);xlim([-1.5 1.5]);ylim([-1.5 1.5]);
grid off

CH1 = convhull(PPim(:,1),PPim(:,2),PPim(:,3));

% Convex Hull과 점들을 플로팅
figure;
trisurf(CH1, PPim(:,1),PPim(:,2),PPim(:,3), 'Facecolor', 'yellow', 'FaceAlpha', 0.95); % Convex Hull의 삼각면을 플로팅
hold on;

axis equal
zlim([-1.5 1.5]);xlim([-1.5 1.5]);ylim([-1.5 1.5]);
hold off
grid off
%}


%%
syms x;
syms y;
syms z;

%% 6각기둥 만들고 그리기
figure
scatter3(PPim(:,1),PPim(:,2),PPim(:,3),3,PPim(:,3),'filled');
xlabel ('X (m)','FontSize',16)
ylabel ('Y (m)','FontSize',16)
zlabel ('Z (m)','FontSize',16)
view([1, 1, 1]);

axis equal
axis([-1 1 -1 1 -1.2 1.2]);

saveas(gcf,'image_ksas\Cuboid_pointcloud_3.svg')
savefig(gcf,'image_ksas\Cuboid_pointcloud_3.fig')

f_cuboid = @(x,y,z)  16*x^4+16*y^4+z^4-1;
figure
fimplicit3(f_cuboid,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
%saveas(gcf,'\image_ksas\Cuboid_poly.svg')

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);

axis equal
axis([-1 1 -1 1 -1.2 1.2]);


saveas(gcf,'image_ksas\Cuboid_fit_3.svg')
savefig(gcf,'image_ksas\Cuboid_fit_3.fig')


%%
f_ellipsoid = @(x,y,z)  2.24*x^2+2.24*y^2+z^2-1;
figure
fimplicit3(f_ellipsoid,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);

axis equal
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

saveas(gcf,'image_ksas\Ellipsoid_fit.svg')
savefig(gcf,'image_ksas\Ellipsoid_fit.fig')

%f_cylinder = @(x,y,z)  (2.24*x^2+2.24*y^2)^2-0.5*(2.24*x^2+2.24*y^2)*z^2+z^4-1;
f_cylinder = @(x,y,z)  (2.24*x^2+2.24*y^2)^2+z^4-1;
figure
fimplicit3(f_cylinder,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);

axis equal
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

saveas(gcf,'image_ksas\Cylinder.svg')
savefig(gcf,'image_ksas\Cylinder.fig')

%%
figure

f_triangle = @(x,y,z)  (x^2+y^2)^3+0.8*(x^3-3*x*y^2)*(x^2+y^2)^1.5-0.5*((x^2+y^2)^1.5+0.8*(x^3-3*x*y^2))*z^3+z^6-1;
%f_triangle = @(x,y,z)  (x^2+y^2)^1.5+0.2*(x^3-3*x*y^2)+z^6-1;
fimplicit3(f_triangle,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);

axis equal
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

saveas(gcf,'image_ksas\Triangle.svg')
savefig(gcf,'image_ksas\Triangle.fig')




%% ngd fitting

figure

f2_ngd = @(x,y,z)  5*( x^4*(y^2+z^2) +y^4*(z^2+x^2)+z^4*(x^2+y^2))-1;
f2_ngd = @(x,y,z)  x^2*y^2*z^2-1;
f2_ngd = @(x,y,z)  x^6+y^6+z^6 - 0.40*( x^4*(y^2+z^2) +y^4*(z^2+x^2)+z^4*(x^2+y^2))+0.9*x^2*y^2*z^2-1;
fimplicit3(f2_ngd,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
saveas(gcf,'\image_ksas\Cuboid_extended2.svg')
f2_ngd = @(x,y,z)  x^4+y^4+z^4-0.5*x^2*z^2-0.5*y^2*z^2-0.5*y^2*x^2-1;
f2_ngd = @(x,y,z)  1.1*x^4+1.1*y^4-1.0*y^2*x^2+z^4-1;
f2_ngd = @(x,y,z)  x^4+y^4-0.0*y^2*x^2-0.0*x^2*z^2-0.0*y^2*z^2+z^4-1;

f2_ngd = @(x,y,z)  x^6 - 15*x^4*y^2+15*x^2*y^4 - y^6-1;
f2_ngd = @(x,y,z)  1.35*(x^2+y^2)^3+z^6 - 0.35*( x^6 - 15*x^4*y^2+15*x^2*y^4 - y^6)-1;

fimplicit3(f2_ngd,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
hold on 
%scatter3(PP01(:,1),PP01(:,2),PP01(:,3),1,PP01(:,3),'filled');
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
hold off
saveas(gcf,'image_ksas\basis_1.svg')
savefig(gcf,'image_ksas\Cuboid_basis_05.fig')


%basis
f2_base = @(x,y)  (x^4+y^4-6*y^2*x^2)*(x^2+y^2);
f2_base = @(x,y)  (4*x*y*(x^2-y^2))*(x^2+y^2);
fcontour(f2_base,[-2 2 -2 2],'Fill', 'on','LineColor','k')
xlabel ('X (m)')
ylabel ('Y (m)')
axis equal
xticks([-2:1:2])
yticks([-2:1:2])
saveas(gcf,'image_ksas\basis_6_2.svg')
savefig(gcf,'image_ksas\basis_6_2.fig')

%%
f7 = x^4+z^4+2*x^2*z^2+0.01*(y-1.2)^4-0.1*(x^2+z^2)*(y-1.2)^2- 0.0001;
f7_partial = matlabFunction(f7, 'Vars', [x, y, z]);
f8 = x^4+z^4+2*x^2*z^2+0.01*(y+1.2)^4-0.1*(x^2+z^2)*(y+1.2)^2- 0.0001;
f8_partial = matlabFunction(f8, 'Vars', [x, y, z]);
f9 = (8*x)^4+(2*z)^4-0.62*(8*x)^2*(2*z)^2+0.66^4*(y-3.0)^4-0.61*0.66^2*(y-3.0)^2*((8*x)^2+(2*z)^2)- 1;
f9_partial = matlabFunction(f9, 'Vars', [x, y, z]);
f10 = (8*x)^4+(2*z)^4-0.62*(8*x)^2*(2*z)^2+0.66^4*(y+3.0)^4-0.61*0.66^2*(y+3.0)^2*((8*x)^2+(2*z)^2)- 1;
f10_partial = matlabFunction(f10, 'Vars', [x, y, z]);

figure
hold on
%fimplicit3(f6_partial,[-1.5 1.5 -1.5 1.5 -1.5 1.5],'FaceColor', 'r');
fimplicit3(f7_partial,[-1.5 1.5 -1.5 3.5 -1.5 1.5],'FaceColor', 'b');
fimplicit3(f8_partial,[-1.5 1.5 -3.5 1.5 -1.5 1.5],'FaceColor', 'b');
fimplicit3(f9_partial,[-1.5 1.5 0.5 4.5 -1.5 1.5],'FaceColor', 'm');
fimplicit3(f10_partial,[-1.5 1.5 -4.5 -0.5 -1.5 1.5],'FaceColor', 'c');
hold off
axis equal
zlim([-3.5 3.5]);xlim([-3.5 3.5]);ylim([-4.5 4.5]);




%%


order = 6;

TermsB = homogeneTerm(order);
betaB = sym('beta', [1, length(TermsB)]);
f2 = sum(betaB .* TermsB);
difx = diff(f2, x);
dify = diff(f2, y);
difz = diff(f2, z);

%0.08s
f2_numeric = matlabFunction(f2, 'Vars', {[x, y, z], betaB});
d1_numeric = matlabFunction(difx, 'Vars', {[x, y, z], betaB});
d2_numeric = matlabFunction(dify, 'Vars', {[x, y, z], betaB});
d3_numeric = matlabFunction(difz, 'Vars', {[x, y, z], betaB});


%[beta_value0,error]    = regressionFourthOrder(PP02);
PPm_use = PPim;
[beta_values,error]    = regressionFourthOrder(PPm_use,TermsB);


f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;



% 초기 값 설정
PPm_shift = PPm_use;  % 초기 PPm 설정
error_shift = error;  % 초기 에러 설정
shiftSum = [0; 0; 0];
shiftResidue = shiftReal.';
shiftResidueSet = [];
numIterations = 7;  % 반복 횟수 설정


figure

fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5],'FaceColor', [0.93 0.83 0.32]);
hold on

scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),2,[0 0.42 0.72],'filled');
hold off
axis equal

view([1 1 1]);
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')



%%
%not useful in this sInario
for i = 1:2*numIterations
    % 미분함수에 정의
    dx = d1_numeric(PPm_shift, beta_values.');
    dy = d2_numeric(PPm_shift, beta_values.');
    dz = d3_numeric(PPm_shift, beta_values.');
    dxyz = [dx dy dz];


    [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);    
    shiftSum = shiftSum-shift;
    shiftResidue = shiftResidue + shift;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    PPm_shift = PPm_shift + shift.';

    [beta_values, error_shift] = regressionFourthOrder(PPm_shift,TermsB);

end

f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
figure
fimplicit3(f2_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5],'FaceColor', [0.93 0.83 0.32]);
hold on

scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),2,[0 0.42 0.72],'filled');
hold off
axis equal

view([1 1 1]);
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')


figure(17)
subplot(1,3,1)
plot(shiftResidueSet(1,:))
subplot(1,3,2)
plot(shiftResidueSet(2,:))
subplot(1,3,3)
plot(shiftResidueSet(3,:))

%%
saveas(gcf,'image_ksas\Iter0.svg')
savefig(gcf,'image_ksas\Iter0.fig')
