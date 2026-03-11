PP05 = generateRandomPointsOnIcosahedron(4000)+randn(4000,3)*0.003;  
PP03 = generateRandomPointsOnCylinder(600)+randn(600,3)*0.003;


q = [0.9239, 0.3827, 0, 0];
Rotm = quat2rotm(q);

PP_cvd = PP03;

PP_cvd(:,1) = PP_cvd(:,1) * 1.276;
PP_cvd(:,2) = PP_cvd(:,2) * 0.8;
PP_cvd(:,3) = PP_cvd(:,3) * 1.65;

PPr_cvd = PP03 * Rotm.';  

cov_PP = cov(PP_cvd);
cov_PPr = cov(PPr_cvd);

[vec_PP, val_PP] = eig(cov_PP);    % PP01의 PCA
[vec_PPr, val_PPr] = eig(cov_PPr); % PPr01의 PCA

figure
center = mean(PP_cvd);
quiver3(center(1), center(2), center(3), vec_PP(1,1), vec_PP(2,1), vec_PP(3,1), 'r')
hold on
quiver3(center(1), center(2), center(3), vec_PP(1,2), vec_PP(2,2), vec_PP(3,2), 'g')
quiver3(center(1), center(2), center(3), vec_PP(1,3), vec_PP(2,3), vec_PP(3,3), 'b')

center = mean(PPr_cvd);
quiver3(center(1), center(2), center(3), vec_PPr(1,1), vec_PPr(2,1), vec_PPr(3,1), 'r')
quiver3(center(1), center(2), center(3), vec_PPr(1,2), vec_PPr(2,2), vec_PPr(3,2), 'g')
quiver3(center(1), center(2), center(3), vec_PPr(1,3), vec_PPr(2,3), vec_PPr(3,3), 'b')
scatter3(PPr_cvd(:,1), PPr_cvd(:,2), PPr_cvd(:,3), '.')

title('PP PCA')
xlabel('x','FontSize',16); ylabel('y','FontSize',16);zlabel('z','FontSize',16);
hold off
axis equal

%%
mean_P = mean(PP_cvd, 1);
mean_Q = mean(PPr_cvd, 1);

X = PP_cvd - mean_P;
Y = PPr_cvd - mean_Q;

cov_cross = (X' * Y)/(length(PP_cvd)-1);
[vec1_PPcr, val_PPcr,vec2_PPcr] = svd(cov_cross);

figure
center = mean(PPr_cvd);
quiver3(center(1), center(2), center(3), vec1_PPcr(1,1), vec1_PPcr(2,1), vec1_PPcr(3,1), 'r')
hold on
quiver3(center(1), center(2), center(3), vec1_PPcr(1,2), vec1_PPcr(2,2), vec1_PPcr(3,2), 'g')
quiver3(center(1), center(2), center(3), vec1_PPcr(1,3), vec1_PPcr(2,3), vec1_PPcr(3,3), 'b')
quiver3(center(1), center(2), center(3), vec2_PPcr(1,1), vec2_PPcr(2,1), vec2_PPcr(3,1), 'r')
hold on
quiver3(center(1), center(2), center(3), vec2_PPcr(1,2), vec2_PPcr(2,2), vec2_PPcr(3,2), 'g')
quiver3(center(1), center(2), center(3), vec2_PPcr(1,3), vec2_PPcr(2,3), vec2_PPcr(3,3), 'b')
%scatter3(PPr_cvd(:,1), PPr_cvd(:,2), PPr_cvd(:,3), '.')
xlabel('x','FontSize',16); ylabel('y','FontSize',16);zlabel('z','FontSize',16);

title('PPcr PCA')
axis equal

%% 다항식의 cvd

syms x y z

f1 =  x^2 + 1.4*y^2 + 2.25*z^2

q = [0.866; 0; 0.5; 0];
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



%% 벡터 랜덤화 검증

quaterinit=zeros(4000, 4);
for i = 1:4000
    qRand = randn(1, 4); 
    qRand = qRand / norm(qRand);
    
    quaterinit(i, :) = qRand;
end

%q1 = [q0; q1; q2; q3];
M = @(q) [ ...
    q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, ...        % x'
    2*(q(2)*q(3) + q(1)*q(4)), ...               % y'
    2*(q(2)*q(4) - q(1)*q(3)) ];

rotated = zeros(size(quaterinit,1), 3);
for i = 1:size(quaterinit,1)
    rotated(i,:) = M(quaterinit(i,:));
end
figure;
scatter3(rotated(:,1), rotated(:,2), rotated(:,3), 10, '.', 'MarkerEdgeAlpha',0.6);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('(1,0,0) rotated by 4000 random quaternions');
grid on;

figure;
histogram(rotated(:,1),'BinLimits', [-1 1],NumBins=20);
%% 랜덤 벡터 and quaternion

mean_G = mean(g_vec, 2);
mean_Q = mean(q_vec, 2);

X = g_vec - mean_G;
Y = q_vec - mean_Q;
cov_X = cov(X.');
[vec_X, val_X] = eig(cov_X);
cov_Y = cov(Y.');
[vec_Y, val_Y] = eig(cov_Y);

figure
quiver3(mean_G(1), mean_G(2), mean_G(3), vec_X(1,1), vec_X(2,1), vec_X(3,1), 'r')
hold on
quiver3(mean_G(1), mean_G(2), mean_G(3), vec_X(1,2), vec_X(2,2), vec_X(3,2), 'g')
quiver3(mean_G(1), mean_G(2), mean_G(3), vec_X(1,3), vec_X(2,3), vec_X(3,3), 'b')
hold off
title('G PCA')
xlabel('x','FontSize',16); ylabel('y','FontSize',16);zlabel('z','FontSize',16);
hold off
axis equal

figure
quiver3(mean_Q(2), mean_Q(3), mean_Q(4), vec_Y(2,1), vec_Y(3,1), vec_Y(4,1), 'r')
hold on
quiver3(mean_Q(2), mean_Q(3), mean_Q(4), vec_Y(2,2), vec_Y(3,2), vec_Y(4,2), 'g')
quiver3(mean_Q(2), mean_Q(3), mean_Q(4), vec_Y(2,3), vec_Y(3,3), vec_Y(4,3), 'b')
quiver3(mean_Q(2), mean_Q(3), mean_Q(4), vec_Y(2,4), vec_Y(3,4), vec_Y(4,4), 'k')
scatter3(q_vec(2,:),q_vec(3,:),q_vec(4,:),1,q_vec(1,:),'.');
hold off

title('Q PCA')
xlabel('x','FontSize',16); ylabel('y','FontSize',16);zlabel('z','FontSize',16);
hold off
axis equal
%%
location = zeros(5,3);
location(1,:) = [-1, 1, 0] ;
location(2,:) = [-1, 0, 1];
location(3,:) = [-1,-1,-1];
location(4,:) = [ 1,-1, 0];
location(5,:) = [ 0, 0, 0];

cov_loc = cov(location);
[vec_loc, val_loc] = eig(cov_loc);

figure(1)
center = [0,0,0];
quiver3(center(1), center(2), center(3), vec_loc(1,1), vec_loc(2,1), vec_loc(3,1), 'r')
hold on
quiver3(center(1), center(2), center(3), vec_loc(1,2), vec_loc(2,2), vec_loc(3,2), 'g')
quiver3(center(1), center(2), center(3), vec_loc(1,3), vec_loc(2,3), vec_loc(3,3), 'b')

title('PP PCA')
xlabel('x','FontSize',16); ylabel('y','FontSize',16);zlabel('z','FontSize',16);
hold off
axis equal