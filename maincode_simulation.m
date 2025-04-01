%%
n  = sqrt(3.986*10^5/(6380+650)^3); 

D1_new = 0;
C2_new = -30;
C3_new = 0;
C4_new = -0;

x0_new = D1_new + C2_new;
y0_new = 2*C3_new+C4_new;
vx0_new = C3_new*n;
vy0_new = -3/2*n*D1_new-2*C2_new;


% (2) 시간 범위 설정 
t = linspace(0, 5400, 100);
x1 = D1_new + C2_new*cos(n*t) + C3_new*sin(n*t);
y1 = -3/2*n*D1_new.*t - 2*C2_new*sin(n*t) + 2*C3_new*cos(n*t) + C4_new;
Orbit = [t(:), x1(:), y1(:),zeros(100,1)];

%
x_max = -20;
x_min = -40;

D1_old = (x_max+x_min)/2;
C_norm = (x_max-x_min)/2;

x0_old = x0_new;
y0_old = y0_new;

C2_old = x0_old - D1_old;
C3_old = sqrt(C_norm^2-C2_old^2);
C4_old = y0_old - 2*C3_old;

t = linspace(-11200, 0, 200);
x1 = D1_old + C2_old*cos(n*t) + C3_old*sin(n*t);
y1 = -3/2*n*D1_old.*t - 2*C2_old*sin(n*t) + 2*C3_old*cos(n*t) + C4_old;
Orbit = [t(:), x1(:), y1(:),zeros(200,1);Orbit];

figure
hold on;
plot(Orbit(:,3), Orbit(:,2), 'LineWidth', 1.2);
set(gca,'XDir','reverse')

hold off;
axis equal;


%%
PP01 = generateRandomPointsOnHexagonPrism(100)+randn(100,3)*0.02;
PP02 = generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 = generateRandomPointsOnCylinder(7000)+randn(7000,3)*0.0005;
PP04 = generateRandomPointsOnCube(200);

PPm_body = PP01 ;
PPm_body(:,1) = PPm_body(:,1) * 0.576;
PPm_body(:,2) = PPm_body(:,2) * 0.576;
PPm_body(:,3) = PPm_body(:,3) * 1.165;

order = 6;
n = [0 sin(35/57.92) cos(35/57.92)];
w = 0.006;
t = linspace(-11200, 0, 200).';
quater0= [cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];

t = linspace(0, 5400, 100).';
quater0= [quater0 ;cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];

%scatter3(quater0(:,2),quater0(:,3),quater0(:,4))


%%

Qinit = [1;0;0;0];
Array_Disp=[];
Array_QB=[];
Array_BetaB=[];

for i= 170:300
    PP01 = generateRandomPointsOnHexagonPrism(100)+randn(100,3)*0.02;
    PPm_body = PP01 ;
    PPm_body(:,1) = PPm_body(:,1) * 0.576;
    PPm_body(:,2) = PPm_body(:,2) * 0.576;
    PPm_body(:,3) = PPm_body(:,3) * 1.165;

    rotm= quat2rotm(quater0(i,:));
    PPmR_body = PPm_body * rotm.';
    PP_use = PPmR_body - [Orbit(i,2), Orbit(i,3),0];
    disp(i)
    tic
    [DispB,QB,BetaB] = PoliNavigationSolver(0,PP_use,6,Qinit);
    toc
    
    Array_Disp=[Array_Disp ; Orbit(i,1), DispB];
    Array_QB=[Array_QB ;  Orbit(i,1), QB.'];
    Array_BetaB=[ Array_BetaB;  Orbit(i,1), BetaB];
    Qinit =  QB;
end

%%
%save('ResultNavigation5.mat',"Array_BetaB","Array_Disp","Array_QB"); 

figure 
plot(Array_Disp(1:5:end,1),Array_Disp(1:5:end,3),'ro','Markersize',3)
hold on 
plot(Orbit(:,1),-Orbit(:,3),'black','LineWidth',1)
hold off
xlabel('Time','FontSize',16); ylabel('Y (m)','FontSize',16);
legend('Estimation','Ground Truth');
saveas(gcf,'image_ksas\Simulation_y.svg')
savefig(gcf,'image_ksas\Simulation_y.fig')
figure 
plot(Array_Disp(1:3:end,1),Array_Disp(1:3:end,2),'ro','Markersize',3)
hold on 
plot(Orbit(:,1),-Orbit(:,2),'black','LineWidth',1)
hold off
xlabel('Time','FontSize',16); ylabel('X (m)','FontSize',16);
legend('Estimation','Ground Truth');
saveas(gcf,'image_ksas\Simulation_x.svg')
savefig(gcf,'image_ksas\Simulation_x.fig')
figure 
plot(Array_Disp(1:3:end,2),Array_Disp(1:3:end,3),'ro','Markersize',3)
hold on 
plot(-Orbit(:,2),-Orbit(:,3),'black','LineWidth',1)
hold off
xlabel('X (m)','FontSize',16); ylabel('Y (m)','FontSize',16);
axis equal
legend('Estimation','Ground Truth');
saveas(gcf,'image_ksas\Simulation_xy.svg')
savefig(gcf,'image_ksas\Simulation_xy.fig')
figure 
plot(Array_Disp(1:131,2)+Orbit(170:300,2),Array_Disp(1:131,3)+Orbit(170:300,3),'ro','Markersize',3)


xlabel('X (m)','FontSize',16); ylabel('Y (m)','FontSize',16);

saveas(gcf,'image_ksas\Simulation_error.svg')
savefig(gcf,'image_ksas\Simulation_error.fig')



figure
scatter3(Array_QB(1:131,3),Array_QB(1:131,4),Array_QB(1:131,5),3,'ro','filled');
hold on
scatter3(quater0(170:300,2),quater0(170:300,3),quater0(170:300,4),3,'k','filled');
hold off
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis equal
legend('Ground Truth','Estimation');

qTemp = zeros(150,4);
for i= 1:131
    qTemp(i,:) = quatmultiply(quatconj(Array_QB(i,2:5)),quater0(i+169,:));
end
scatter3(qTemp(:,2),qTemp(:,3),qTemp(:,4),4,qTemp(:,1),'filled')
axis equal
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis ([-0.5 0.5 -0.5 0.5 -1 1]);
for i= 1:50
    qTemp2(i,:) = quatmultiply(quatconj(qTemp(i,1:4)),qTemp(6,1:4));
end
figure
scatter3(qTemp2(:,2),qTemp2(:,3),qTemp2(:,4),4,qTemp2(:,1),'filled')
axis equal
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);

%%
TermsB = homogeneTerm(order);
%fA = BetaA *  TermsB(:);
fB = BetaB *  TermsB(:);
%fA_shift =  subs(fA, [x, y, z], [x - DispA(1), y - DispA(2), z - DispA(3)])-1;
fB_shift =  subs(fB, [x, y, z], [x - DispB(1), y - DispB(2), z - DispB(3)])-1;

%{
figure(4)
scatter3(PP_use(:,1),PP_use(:,2),PP_use(:,3),1,PP_use(:,3),'filled');
hold on
fimplicit3(fB_shift,[-100 100 -1000 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
grid on
%}
%%
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
fB_rot = subs(fB_rot, [q0,q1,q2,q3], QB.');
[coeffsfB, monomialfB] = coeffs(fB_rot , [x,y,z]);
coeffsfB = double(coeffsfB);
coeffsfB_unused = sym([]);
monList = [
  0 6 0;
  2 4 0;   
  4 2 0;
  6 0 0;   
  0 0 6;
];

for k = 1:length(monomialfB) %안쓰는 항들의 계수 norm 구하기
    degx = feval(symengine, 'degree', monomialfB(k), x);
    degy = feval(symengine, 'degree', monomialfB(k), y);
    degz = feval(symengine, 'degree', monomialfB(k), z);
    if ~ismember([degx, degy, degz], monList, 'rows');
       coeffsfB_unused(end+1) = coeffsfB(k); 
       [degx, degy, degz]
    end
end



%RqA = double(subs(Rq,[q0,q1,q2,q3],[QA(1),-QA(2),-QA(3),-QA(4)]));
%PPmR_shiftA = PPmR_body - DispA;
%PPmR_rotA = (RqA * PPmR_shiftA.').';

RqB = double(subs(Rq,[q0,q1,q2,q3],[QB(1),-QB(2),-QB(3),-QB(4)]));
PPmR_shiftB = PPmR_body;
PPmR_rotB = (RqB * PPmR_shiftB.').';
%{
figure(1)
scatter3(PPmR_rotA(:,1),PPmR_rotA(:,2),PPmR_rotA(:,3),1,PPmR_rotA(:,3),'filled');
hold on
fimplicit3(fA_rot-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid on
%}
figure(1)
scatter3(PPmR_body(:,1),PPmR_body(:,2),PPmR_body(:,3),3,PPmR_body(:,3),'filled');
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal

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
figure(2)
scatter3(PPmR_rotB(:,1),PPmR_rotB(:,2),PPmR_rotB(:,3),3,PPmR_rotB(:,3),'filled');
hold on
fimplicit3(fB_rot-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid on

