function output = calculateDerivative(inputMatrix)
    % 입력 행렬에서 각 행을 추출하여 각각의 변수 x, y, z로 지정
    x = inputMatrix(:, 1);
    y = inputMatrix(:, 2);
    z = inputMatrix(:, 3);
    
    % 각 항을 계산
    dx = 4*x.^3;
    dy = 4*y.^3;
    dz = 4*z.^3;
    
    
    % 결과 행렬을 생성
    output = [dx dy dz];
end