%symbolic test

syms x y z
beta = sym('beta', [1, 15]);  % 15개의 심볼릭 변수로 beta 정의
tic
% 예시로 f(x, y, z) 함수 정의 (심볼릭 beta 사용)
f = beta(1)*x^4 + beta(2)*y^4 + beta(3)*z^4 + ...
    beta(4)*x^2*y^2 + beta(5)*x^2*z^2 + beta(6)*y^2*z^2 + ...
    beta(7)*x^3*y + beta(8)*x^3*z + beta(9)*y^3*x + ...
    beta(10)*y^3*z + beta(11)*z^3*x + beta(12)*z^3*y - 1 + ...
    beta(13)*x^2*y*z + beta(14)*y^2*z*x + beta(15)*z^2*x*y;
toc
% 함수 재정의에 0.03초 소요
%%
% f에 대해 x, y, z에 대한 미분 수행

difx = diff(f, x);
dify = diff(f, y);
difz = diff(f, z);

%%
% 미분식을 한번만 계산한 후, 수치적 계산에 재사용
tic
difx_func = matlabFunction(difx, 'Vars', {[x, y, z], beta});
dify_func = matlabFunction(dify, 'Vars', {[x, y, z], beta});
difz_func = matlabFunction(difz, 'Vars', {[x, y, z], beta});
toc
%병목 함수 변환에 0.2초 소요. 
%%
% 이후 수치적으로 계산할 때는 beta 값을 바꿔가며 사용
% 예시로 beta 값을 수치적으로 정의

for i = 1:5
    beta_values = rand(1, 15);  % 15개의 수치적 값으로 beta 정의
    
    
    
    % 좌표 예시 (1, 2, 3)
    coords = [1, 1/sqrt(2), 1/sqrt(2)];
    tic
    % 미분값 계산
    dx = difx_func(PP, beta_values);
    dy = dify_func(PP, beta_values);
    dz = difz_func(PP, beta_values);
    toc
    % 결과 출력
end
% 대입 한번에  0.002초 1000번 하면 2초