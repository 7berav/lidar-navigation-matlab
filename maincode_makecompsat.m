%2024 08 20
PP01 = generateRandomPointsOnHexagonPrism(30000)+randn(30000,3)*0.0007;
PP02 = generateRandomPointsOnCube(1000)+randn(1000,3)*0.0005;
PP03 = generateRandomPointsOnCylinder(7000)+randn(7000,3)*0.0005;
PP04 = generateRandomPointsOnCube(200);
shiftReal =  [-3 -0 56];
PPm_body = PP01 ;
PPm_body(:,1) = PPm_body(:,1) * 0.576;
PPm_body(:,2) = PPm_body(:,2) * 0.576;
PPm_body(:,3) = PPm_body(:,3) * 1.165;
PPm_body = PPm_body(PPm_body(:,3)-1.5*PPm_body(:,1)>-0.6,:); 
PPm_body = PPm_body + shiftReal;

PPm_cyl = PP03 ;
PPm_cyl(:,1) = PPm_cyl(:,1) * 0.6;
PPm_cyl(:,2) = PPm_cyl(:,2) * 0.6;
PPm_cyl(:,3) = PPm_cyl(:,3) * 1.165;
PPm_cyl = PPm_cyl(PPm_cyl(:,2)>-0.05&PPm_cyl(:,3)>-1.1,:); 
PPm_cyl = PPm_cyl + shiftReal;
%Vector_temp = randn(size(PPm_body,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.3;
%PPm_body = PPm_body + Vector_temp;

%{
PPm_panel1 = PP02 ;
PPm_panel1(:,1) = PPm_panel1(:,1) / 16;
PPm_panel1(:,2) = PPm_panel1(:,2) * 1.25;
PPm_panel1(:,3) = PPm_panel1(:,3) * 0.576;
PPm_panel1 = PPm_panel1 + [0 2.025 0];
%Vector_temp = randn(size(PPm_panel1,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_panel1 = PPm_panel1 + Vector_temp;


PPm_panel2 = PP02 ;
PPm_panel2(:,1) = PPm_panel2(:,1) / 16;
PPm_panel2(:,2) = PPm_panel2(:,2) * 1.25;
PPm_panel2(:,3) = PPm_panel2(:,3) * 0.576;
PPm_panel2 = PPm_panel2 - [0 2.025 0];
%Vector_temp = randn(size(PPm_panel2,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_panel2 = PPm_panel2 + Vector_temp;

PPm_box1 = PP04 ;
PPm_box1(:,1) = PPm_box1(:,1) * 0.333;
PPm_box1(:,2) = PPm_box1(:,2) * 0.175;
PPm_box1(:,3) = PPm_box1(:,3) * 0.2;
PPm_box1 = PPm_box1 + [0 +0.295 -1.365];
%Vector_temp = randn(size(PPm_box1,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_box1 = PPm_box1 + Vector_temp;

PPm_EOC = PP03 ;
PPm_EOC(:,1) = PPm_EOC(:,1) * 0.09;
PPm_EOC(:,2) = PPm_EOC(:,2) * 0.09;
PPm_EOC(:,3) = PPm_EOC(:,3) * 0.05;
PPm_EOC = PPm_EOC + [-0.115 +0.295 -1.615];
%Vector_temp = randn(size(PPm_EOC,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_EOC = PPm_EOC + Vector_temp;

PPm_box2 = PP04 ;
PPm_box2(:,1) = PPm_box2(:,1) * 0.05;
PPm_box2(:,2) = PPm_box2(:,2) * 0.2;
PPm_box2(:,3) = PPm_box2(:,3) * 0.1;
PPm_box2 = PPm_box2 + [0.2 -0.215 -1.265];
%Vector_temp = randn(size(PPm_box2,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_box2 = PPm_box2 + Vector_temp;

PPm_box3 = PP04 ;
PPm_box3(:,1) = PPm_box3(:,1) * 0.05;
PPm_box3(:,2) = PPm_box3(:,2) * 0.2;
PPm_box3(:,3) = PPm_box3(:,3) * 0.2;
PPm_box3 = PPm_box3 + [-0.2 -0.215 -1.365];
%Vector_temp = randn(size(PPm_box3,1),3)*1;
%Vector_temp = Vector_temp ./ vecnorm(Vector_temp,2,2) * 0.2;
%PPm_box3 = PPm_box3 + Vector_temp;

%PPmm = PPmm + shiftReal;
%PPm_total=[PPm_body ; PPm_panel1; PPm_panel2;PPm_box1;PPm_EOC;PPm_box2;PPm_box3];
PPm_total=[PPm_body ; PPm_box1;PPm_EOC;PPm_box2;PPm_box3];
%}
%q= [1 0 0 0];
qo= [cos(17/57.92) sin(17/57.92)*cos(45/57.92) 0 sin(17/57.92)*sin(45/57.92)];
rotm= quat2rotm(qo);
PPmR_body = PPm_body * rotm.';


%{
figure
h = PPm_total(:,3);
scatter3(PPm_total(:,1),PPm_total(:,2),PPm_total(:,3),1,h,'filled');
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
view([1,1,1])
xlim([-2 2]);  
ylim([-4 4]);
zlim([-2 2]);
%title('Pointcloud of KOMPSAT-1 Model')
%}

%% 
syms x;
syms y;
syms z;
order = 6;
PPm_use = PPmR_body;
varName='body';
%% 그냥 좌표 평균 찾기 겸 점 초기세팅

center_shift = mean(PPm_use,1)
PPm_shift = PPm_use - center_shift;

shiftResidue = shiftReal.';
shiftResidue = shiftResidue - center_shift.';
%% 비동차항 먼저


TermsA = nonhomogeneTerm(order);
betaA = sym('beta', [1, length(TermsA)]);
f1 = sum(betaA .* TermsA);

%점 입력

[beta_values0,error0]    = regressionFourthOrder(PPm_use,TermsA);
f1_total   = subs(f1,betaA,beta_values0.');
syms a b c 
f1_shift = subs(f1_total, [x, y, z], [x+a, y+b, z+c]);
f1_expanded = expand(f1_shift);
[coeffsA, monomialA] = coeffs(f1_shift , [x,y,z]);

coeffsA_unused = sym([]);
for k = 1:length(monomialA)
    deg = feval(symengine, 'degree', monomialA(k), x)+feval(symengine, 'degree', monomialA(k), y)+feval(symengine, 'degree', monomialA(k), z);
    if deg ~= order & deg ~= 0
       coeffsA_unused(end+1) = coeffsA(k);  
    end
end

f_coeffsA_norm = sum(coeffsA_unused.^2);  % f_obj(a,b,c)



%%
initGuess = [0,0,0];
gradF = gradient(f_coeffsA_norm, [a, b, c]);  
fNum = matlabFunction(gradF, 'Vars', [a, b, c]);
fHandle2 = @(var) fNum(var(1), var(2), var(3));
optionsB = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 3000);
tic;
[xSol_A, fval_A, exitflag_B, output_B] = fsolve(fHandle2, initGuess, optionsB);
toc;

minVal = double(subs(f_coeffsA_norm, [a,b,c], xSol_A));
fprintf('해당 해에서의 2-norm^2(비동차항) = %.6f\n', minVal);


PPm_shiftA = PPm_use - xSol_A;
f1_substituted = subs(f1_expanded, [a, b, c], xSol_A);
f1_numeric = matlabFunction(f1_substituted-1, 'Vars', [x, y, z]);

f1_origin_numeric = matlabFunction(f1_total-1, 'Vars', [x, y, z]);
[coeffsA2, monomialA2] =coeffs(f1_substituted-1);
coeffsA2=double(coeffsA2);
%{
figure


%h = PPm_shiftA(:,3);

%fimplicit3(f1_numeric,[-4.5 4.5 -4.5 4.5 -4.5 54.5],'FaceColor', [0.99 0.75 0.12]);
%scatter3(PPm_shiftA(:,1),PPm_shiftA(:,2),PPm_shiftA(:,3),1,'h','filled');
scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),1,PPm_use(:,3),'filled');
hold on
fimplicit3(f1_origin_numeric,[-6 2 -3.5 3.5 53.5 59.5],'FaceColor', [0.99 0.75 0.12]);

hold off
colormap(jet);
view([1,1,1])
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
axis equal
grid on
%}

%% Iteration 초기 설정

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

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

shift_1= [0;0;0];
shiftSet = [];
shift1Set = [];
shiftResidueSet = [0;0;0];
shiftSum = [0 ; 0 ; 0];
numIterations = 20;  % 반복 횟수 설정
shiftSum = center_shift.';




%% 그림 그리기_initial  
f2_partial = @(x,y,z) f2_numeric([x,y,z], beta_values.')-1;
f2_substituted = subs(f2, betaB, beta_values.');
f2_translated = subs(f2_substituted-1, [x,y,z], [x-center_shift(1), y-center_shift(2), z-center_shift(3)]);%여기서 1 뺌
f2_translated_expanded = expand(f2_translated);
f2_handle = matlabFunction(f2_translated_expanded, 'Vars', [x,y,z]);
[coeffs1, monomial] = coeffs(f2_partial , [x,y,z]);
%coeffs(f2_substituted , [x,y,z])

%{
figure
h2 = scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),1,PPm_use(:,3),'filled');
colormap(jet);
hold on
h1 = fimplicit3(f2_translated_expanded,[-6 6 -6 6 -6 6],'FaceColor', [0.99 0.75 0.12]);
hold off

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
%title('Pointcloud Model - Panel')
%}





%%
for i = 1:numIterations
    % 미분함수에 정의
    dx = d1_numeric(PPm_shift, beta_values.');
    dy = d2_numeric(PPm_shift, beta_values.');
    dz = d3_numeric(PPm_shift, beta_values.');
    dxyz = [dx dy dz];

    [shift, residual] = regressionShift(PPm_shift, error_shift,dxyz);% 함수 결과값부터 부호가 반대
    shift_1 = shift*0.85 + shift_1*0.66;
    shift1Set= [shift1Set shift_1];    
    shiftSet = [shiftSet shift];
    shiftResidue = shiftResidue + shift_1;% 부호가 반대인데 알아서 해석하세요
    shiftSum = shiftSum - shift_1;
    shiftResidueSet = [shiftResidueSet shiftResidue];
    
    PPm_shift = PPm_shift + shift_1.';
    [beta_values, error_shift] = regressionFourthOrder(PPm_shift,TermsB);

    %0.001s 
    
end

f3_substituted = subs(f2, betaB, beta_values.');
f3_translated = subs(f3_substituted-1, [x,y,z], [x-shiftSum(1), y-shiftSum(2), z-shiftSum(3)]);
f3_translated_expanded =expand(f3_translated);

%{
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
%}

figure
h2 = scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),1,PPm_shift(:,3),'filled');
colormap(jet);
hold on
h1 = fimplicit3(f3_substituted-1,[-6 6 -6 6 -6 6],'FaceColor', [0.05, 0.2, 0.5]);
hold off

xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);

axis equal
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])

equation2string(beta_values, TermsB);
sum(beta_values);

f2_handle(0,+0.866,0)+1;
f2_handle(1,0,0)+1;
f2_handle(-0.5,-0.866,0)+1;
f2_handle(-0.5,-0.866,1)+1;
%norm(error_shift,1)
%evaluateModel(PPm_shift,f2_numeric, beta_values)
%

%saveas(gcf,'image_ksas\RotHexagon_fit.svg')
%savefig(gcf,'image_ksas\RotHexagon_fit_.fig')

figure(1)
subplot(1,3,1)
plot(shiftResidueSet(1,:))
subplot(1,3,2)
plot(shiftResidueSet(2,:))
subplot(1,3,3)
plot(shiftResidueSet(3,:))
%% 
syms q0 q1 q2 q3 real

% 간단히 벡터화
q = [q0; q1; q2; q3];

% 회전행렬 R(q) 정의
Rq = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
       2*(q2*q1 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
       2*(q3*q1 - q0*q2),   2*(q3*q2 + q0*q1),     q0^2 - q1^2 - q2^2 + q3^2 ];

xr = Rq(1,1)*x + Rq(1,2)*y + Rq(1,3)*z;
yr = Rq(2,1)*x + Rq(2,2)*y + Rq(2,3)*z;
zr = Rq(3,1)*x + Rq(3,2)*y + Rq(3,3)*z;
tic
f4_rotated = subs(f3_substituted, [x,y,z], [xr,yr,zr]);%중심에 있는 f3를 회전
toc
%tic
%f4_expanded = expand(f4_rotated);
%toc
[coeffsB, monomialB] = coeffs(f4_rotated , [x,y,z]);
coeffsB_unused = sym([]);
monList = [
  0 6 0;
  2 4 0;   
  4 2 0;
  6 0 0;   
  0 0 6;
];

for k = 1:length(monomialB) %안쓰는 항들의 계수 norm 구하기
    degx = feval(symengine, 'degree', monomialB(k), x);
    degy = feval(symengine, 'degree', monomialB(k), y);
    degz = feval(symengine, 'degree', monomialB(k), z);
    if ~ismember([degx, degy, degz], monList, 'rows');
       coeffsB_unused(end+1) = coeffsB(k);  
    end
end

f_coeffsB_norm = sum(coeffsB_unused.^2);  % f_obj(a,b,c)

fNum = matlabFunction(f_coeffsB_norm, 'Vars', [q0, q1, q2, q3]);
fHandle3 = @(var) fNum(var(1), var(2), var(3), var(4));
nonlcon = @(qVec) deal([],qVec'*qVec - 1);  


initQ = [1; 0; 0; 0];
optionsC = optimoptions('fmincon','Display','iter');

[qOpt, fValB] = fmincon(@(qIn) fHandle3 (qIn), ...
                       initQ,[],[],[],[],[],[], nonlcon, optionsC)
qOpt = qOpt / norm(qOpt);  % safety normalize


%대입
f4_substituted = subs(f4_rotated, [q0,q1,q2,q3], qOpt.');%중심에 있는 f3를 회전
Rq4 = double(subs(Rq,[q0,q1,q2,q3],[qOpt(1),-qOpt(2),-qOpt(3),-qOpt(4)]));
PPm_shiftB = (Rq4 * PPm_shift.').'; 
f4_numeric = matlabFunction(f4_substituted-1, 'Vars', [x, y, z]);
[coeffsB2, monomialB2] = coeffs(f4_substituted , [x,y,z]);
coeffsB2 = double(coeffsB2);

figure
hold on
scatter3(PPm_shiftB(:,1),PPm_shiftB(:,2),PPm_shiftB(:,3),1,PPm_shiftB(:,3),'filled');
fimplicit3(f4_numeric,[-2 2 -2 2 -2 2],'FaceColor', [0.05, 0.2, 0.5]);
hold off
colormap(jet);
xlabel ('X (m)')
ylabel ('Y (m)')
zlabel ('Z (m)')
view([1, 1, 1]);
axis equal
grid on
