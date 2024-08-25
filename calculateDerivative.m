function output = calculateDerivative(inputMatrix,f)
    % 입력 행렬에서 각 행을 추출하여 각각의 변수 x, y, z로 지정
    x = inputMatrix(:, 1);
    y = inputMatrix(:, 2);
    z = inputMatrix(:, 3);
    n = length(inputMatrix);
    output = zeros(n,3);
    % 각 항을 계산
    %dx = 4*x.^3;
    %dy = 4*y.^3;
    %dz = 4*z.^3;
    difx = diff(f,x);
    dify = diff(f,y);
    difz = diff(f,z);
    
    for i = 1:n
        output (i,1) = double(subs(difx,[x,y,z],inputMatrix(i,:)));
        output (i,2) = double(subs(dify,[x,y,z],inputMatrix(i,:)));
        output (i,3) = double(subs(difz,[x,y,z],inputMatrix(i,:)));
    end
end