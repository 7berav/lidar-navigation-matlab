function output = calculateDerivative(inputMatrix,f)
    % 입력 행렬에서 각 행을 추출하여 각각의 변수 x, y, z로 지정
    
    xval = inputMatrix(:, 1);
    yval = inputMatrix(:, 2);
    zval = inputMatrix(:, 3);
    n = length(inputMatrix);
    output = zeros(n,3);
    
    syms x y z;

    % 각 항을 계산
    %dx = 4*x.^3;
    %dy = 4*y.^3;
    %dz = 4*z.^3;
    
    difx = matlabFunction(diff(f,x), 'Vars', [x, y, z]);
    dify = matlabFunction(diff(f,y), 'Vars', [x, y, z]);
    difz = matlabFunction(diff(f,z), 'Vars', [x, y, z]);
    

    
    dx  = difx(xval,yval,zval);
    dy  = dify(xval,yval,zval);
    dz  = difz(xval,yval,zval);
    
    output = [dx dy dz];
end