%%
PP01 = generateRandomPointsOnHexagonPrism(2000)+randn(2000,3)*0.02;
PP02 = generateRandomPointsOnCube(150)+randn(150,3)*0.02;
PP03 = generateRandomPointsOnCube(15000)+randn(15000,3)*0.005;

PPm_body = PP01 ;
PPm_body(:,1) = PPm_body(:,1) * 0.6;
PPm_body(:,2) = PPm_body(:,2) * 0.6;
PPm_body(:,3) = PPm_body(:,3) * 1.2;


PP03(:,1) = PP03(:,1) * 0.2;
PP03(:,2) = PP03(:,2) * 0.3;
PP03(:,3) = PP03(:,3) * 0.1;

%가시화 변수
camPos = [-20 0 0]; 
eps_r = 1e-9;

n = [0 sin(25/57.92) cos(25/57.92)];
w = 0.01;
t = linspace(0, 5400, 900).';
dt = 5400/900;
q0= [cos(15/57.92), -sin(15/57.92),0 ,0];
quater0 = [cos(w*t), sin(w*t)*n(1), sin(w*t)*n(2), sin(w*t)*n(3)];
for i= 1:900
    quater1(i,:) = quatmultiply(quater0(i,:),q0);
end
ran = 1:5:100;

%% pointcloud 생성
syms x y z
order = 4;
TermsB = homogeneTerm(order);
FuncsB = matlabFunction(TermsB);
Qinit = [1;0;0;0];
Array_Disp=zeros(length(quater0),4);
Array_QB=zeros(length(quater0),5);
Array_BetaB=zeros(length(quater0),1+length(TermsB));

filename = 'PPmR_body_animation2.gif';
[beta_values,error]  = regressionFourthOrder(PPm_body,FuncsB);
fA= TermsB * beta_values;
figure(10)
scatter3(PPm_body(:,1),PPm_body(:,2),PPm_body(:,3),2,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.5);
hold on
fimplicit3(fA-1,[-100 100 -100 100 -100 100],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor', 'none');
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal   
axis ([-1.5 1.5 -1.5 1.5 -1.5 1.5] )  
lighting gouraud; camlight headlight; material dull;
set(gcf,'Renderer','opengl');       % 투명도는 OpenGL 권장
set(gca,'SortMethod','depth');      % 깊이 기준 정렬
%set(gca,'TwoSidedLighting','on');   % 양면 조명
% 개별 표면(핸들 hSurf가 있을 때)

%hSurf = findobj(gca, 'Type', 'Surface');
%set(hSurf,'BackFaceLighting','reverselit');  % 뒷면도 빛 받게

%%
%{
P = PPm_body;          % N x 3
K = convhulln(P);


viewDir = [-1 0 0];                 % 카메라에서 물체로 향하는 방향
viewDir = viewDir / norm(viewDir);  % 단위벡터로 정규화

v1 = P(K(:,2),:) - P(K(:,1),:);   % 첫 번째 에지벡터 (NF x 3)
v2 = P(K(:,3),:) - P(K(:,1),:);   % 두 번째 에지벡터 (NF x 3)
n  = cross(v1, v2, 2);            % 법선벡터 (NF x 3)
n  = n ./ vecnorm(n, 2, 2);       % 각 행을 단위벡터로 정규화

% 4) 앞면 패싯 판별:
frontFace = (n * viewDir.' < 0);    
visibleIdx = unique(K(frontFace,:));  % 정렬된 유일 index 집합
PP_visible = P(visibleIdx,:);         % (-1,0,0)에서 보이는 hull 정점들
maskVisible = false(size(P,1),1);
maskVisible(visibleIdx) = true;



cols = lines(2);        % cols(1,:), cols(2,:) 가 각각 RGB (1x3)
col_vis = cols(2,:);    % 보이는 점에 쓸 색 (1x3 실수벡터, [R G B])

figure;
scatter3(P(:,1), P(:,2), P(:,3), ...
         5, [0.8 0.8 0.8], ...       % 작은 크기, 연한 회색
         'filled', 'MarkerEdgeColor', 'none');
hold on;
grid on;
axis equal;
scatter3(PP_visible(:,1), PP_visible(:,2), PP_visible(:,3), ...
         25, col_vis, ...            % 더 크게, 팔레트 2번 색
         'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('All points (P)', 'Visible points (PP\_visible)');
hold off;

%%
P = PPm_body;    
camPos = [-100 0 0];               % 관측점
% 점군 스케일에 비례한 반전 반지름 설정
R = 2 * max(pdist(P));             % 점군 대략 크기의 2배

% 1) 관측점 기준 벡터
V = P - camPos;
dist2 = sum(V.^2, 2);
invFactor = (R^2) ./ dist2;

% 2) 구면 반전된 점 좌표
P_inv = camPos + V .* invFactor;

% 3) 반전된 점들의 convex hull 계산
K = convhulln(P_inv);
visibleIdx = unique(K(:));
PP_visible = P(visibleIdx,:);
%}
%%
P = generateRandomPointsOnHexagonPrism(2000);                            % N x 3
camPos1 = [-20 0 3];                   % 관측점 C
U = P - camPos1;                        % 각 점까지의 벡터 (N x 3)
r = sqrt(sum(U.^2, 2));                % 각 점까지의 거리 r_i (N x 1)
R = 1.1 * max(r);                      % 반전 반경 R > max(r_i)

r_safe = max(r, eps_r);
U_hat = U ./ r_safe;                   % N x 3
P_inv = U + 2*(R - r_safe).*U_hat;     % viewpoint를 원점으로 본 좌표계에서의 반전 점
P_inv_all = [P_inv; 0 0 0];            % (N+1) x 3, 마지막 행이 원점
K = convhulln(P_inv_all);

hullIdxAll = unique(K(:));             % hull 위의 vertex 인덱스들 (1..N+1)
N = size(P,1);
visibleIdx = hullIdxAll(hullIdxAll <= N);
P = P + randn(2000,3)*0.02;
PP_visible = P(visibleIdx, :);

cols = lines(2);        % cols(1,:), cols(2,:) 가 각각 RGB (1x3)
col_vis = cols(2,:);    % 보이는 점에 쓸 색 (1x3 실수벡터, [R G B])

figure;
scatter3(P(:,1), P(:,2), P(:,3), ...
         5, [0.8 0.8 0.8], ...       % 작은 크기, 연한 회색
         'filled', 'MarkerEdgeColor', 'none');
hold on;
grid on;
axis equal;
scatter3(PP_visible(:,1), PP_visible(:,2), PP_visible(:,3), ...
         25, col_vis, ...            % 더 크게, 팔레트 2번 색
         'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('All points (P)', 'Visible points (PP\_visible)');
hold off;

%%
Binit = beta_values.';
tic
for i= ran
    PP01 = generateRandomPointsOnHexagonPrism(2000);
    PPm_body = PP01 ;
    PPm_body(:,1) = PPm_body(:,1) * 0.6;
    PPm_body(:,2) = PPm_body(:,2) * 0.6;
    PPm_body(:,3) = PPm_body(:,3) * 1.2;

    rotm= quat2rotm(quater1(i,:));
    PPmR_body = PPm_body * rotm.';
    PP_use = PPmR_body;

    
    P = PP_use;          % N x 3
    U = P - camPos;                        % 각 점까지의 벡터 (N x 3)
    r = sqrt(sum(U.^2, 2));                % 각 점까지의 거리 r_i (N x 1)
    R = 1.1 * max(r);                      % 반전 반경 R > max(r_i)
    
    r_safe = max(r, eps_r);
    U_hat = U ./ r_safe;                   % N x 3
    P_inv = U + 2*(R - r_safe).*U_hat;     % viewpoint를 원점으로 본 좌표계에서의 반전 점
    P_inv_all = [P_inv; 0 0 0];            % (N+1) x 3, 마지막 행이 원점
    K = convhulln(P_inv_all);
    
    hullIdxAll = unique(K(:));             % hull 위의 vertex 인덱스들 (1..N+1)
    N = size(P,1);
    visibleIdx = hullIdxAll(hullIdxAll <= N);
    
    P = P + randn(2000,3)*0.02;
    PP_visible = P(visibleIdx, :);
    
    PP_use = PP_visible;


    disp(i)
    tic
    [DispB,QB,BetaB] = PoliNavigationSolver(0,PP_use,order,Binit,Qinit);
    toc
    
    fA= TermsB * BetaB.';
    figure(10); clf;
    
    scatter3(PPmR_body(:,1), PPmR_body(:,2), PPmR_body(:,3), 4,[0.6 0.6 0.6], 'filled');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    hold on
    scatter3(PP_use(:,1), PP_use(:,2), PP_use(:,3), 4, 'filled');
    fimplicit3(fA-1,[-10 10 -10 10 -10 10],'FaceColor', [0.95, 0.82, 0.5]);
    hold off
    set(gcf, 'Color', 'w');
    grid on;
    xlim([-4 4]); ylim([-4 4]); zlim([-4 4]);  % 범위는 적절히 조정
    %axis ([-4 4 -4 4 -4 4] )  
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    %{
    % GIF에 저장
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
    %}
    Array_Disp(i,:) = [t(i,1), DispB];
    Array_QB(i,:) = [t(i,1), QB.'];
    Array_BetaB(i,:) =[t(i,1), BetaB];
    Qinit =  QB;
end
toc
%% result
%save('ResultQuaternion2.mat',"Array_BetaB","Array_Disp","Array_QB"); 


figure(11)
scatter3(Array_Disp(ran,2),Array_Disp(ran,3),Array_Disp(ran,4),3,'bo','filled')
xlabel('X (m)','FontSize',16); ylabel('Y (m)','FontSize',16);

%saveas(gcf,'image_ksas\Simulation_error.svg')
%savefig(gcf,'image_ksas\Simulation_error.fig')


figure(12)
scatter3(Array_QB(ran,3),Array_QB(ran,4),Array_QB(ran,5),3,'ro','filled');
hold on
scatter3(quater1(ran,2),quater1(ran,3),quater1(ran,4),3,'k','filled');
hold off
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis equal
xlim([-0.5 0.5]); ylim([-0.5 0.5]);
legend('Estimation','Ground Truth');

figure(13)
plot(Array_QB(ran,1),Array_QB(ran,3),'r-');
hold on
plot(Array_QB(ran,1),quater1(ran,2),'k-');

%plot(Array_QB(:,1),Array_QB(:,4),'g-');
%plot(Array_QB(:,1),quater1(range,3),'b-');
%plot(Array_QB(:,1),Array_QB(:,5),'ro');
%plot(Array_QB(:,1),quater1(range,4),'ko');

hold off
xlabel('t','FontSize',16); ylabel('q1','FontSize',16);
legend('Estimation','Ground Truth');


qTemp = zeros(length(quater0),4);
for i= ran
    qTemp(i,:) = quatmultiply(quatconj(Array_QB(i,2:5)),quater1(i,:));
end
figure(14);
scatter3(qTemp(ran,2),qTemp(ran,3),qTemp(ran,4),4,qTemp(ran,1),'filled')
axis equal
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);
axis ([-0.5 0.5 -0.5 0.5 -1 1]);
%%
figure;
histogram(qTemp(ran,4), 20);  % q3 성분
grid on;
xlabel('q_{e3}');
ylabel('Count');
title('Error quaternion q3 histogram');
q0 = qTemp(ran, 1);  
theta = 2*acos( min(1, abs(q0)) );  % [rad]
theta_deg = theta * 180/pi;        % [deg]

figure;
histogram(theta_deg, 30);
grid on;
xlabel('Rotation error \theta (deg)');
ylabel('Count');
title('Rotation error angle histogram');
%%

for i= 1:50
    qTemp2(i,:) = quatmultiply(quatconj(qTemp(i,1:4)),qTemp(6,1:4));
end
figure
scatter3(qTemp2(:,2),qTemp2(:,3),qTemp2(:,4),4,qTemp2(:,1),'filled')
axis equal
xlabel('q1','FontSize',16); ylabel('q2','FontSize',16);zlabel('q3','FontSize',16);

%%
TermsB = homogeneTerm(order);

fB = BetaB *  TermsB(:);

fB_shift =  subs(fB, [x, y, z], [x - DispB(1), y - DispB(2), z - DispB(3)])-1;


%% rotated image making
syms q0 q1 q2 q3 real
syms x y z
Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
       2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
       2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];

xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;



fB_rot = subs(fB, [x,y,z], [xr,yr,zr]);
fB_rot = subs(fB_rot, [q0,q1,q2,q3], QB.');
[coeffsfB, monomialfB] = coeffs(fB_rot , [x,y,z]);
coeffsfB = double(coeffsfB);
coeffsfB_unused = sym([]);



RqB = double(subs(Rq,[q0,q1,q2,q3],[QB(1),-QB(2),-QB(3),-QB(4)]));
PPmR_shiftB = PPmR_body;
PPmR_rotB = (RqB * PPmR_shiftB.').';

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
scatter3(PP03(:,1),PP03(:,2),PP03(:,3),2,'blue');


hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid on

