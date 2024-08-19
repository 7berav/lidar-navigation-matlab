function output = calculateFourthOrder(inputMatrix)
    % 입력 행렬에서 각 행을 추출하여 각각의 변수 x, y, z로 지정
    x = inputMatrix(:, 1);
    y = inputMatrix(:, 2);
    z = inputMatrix(:, 3);
    
    % 각 항을 계산
    x4 = x.^4;
    y4 = y.^4;
    z4 = z.^4;
    x2y2 = (x.^2).*(y.^2);
    x2z2 = (x.^2).*(z.^2);
    y2z2 = (y.^2).*(z.^2);
    x3y = (x.^3).*y;
    x3z = (x.^3).*z;
    y3x = (y.^3).*x;
    y3z = (y.^3).*z;
    z3x = (z.^3).*x;
    z3y = (z.^3).*y;

    x2yz = (x.^2).*y.*z;
    xy2z = (y.^2).*z.*x;
    xyz2 = (z.^2).*x.*y;
    
    
    % 결과 행렬을 생성
    output = [x4 y4 z4 x2y2 x2z2 y2z2 x3y x3z y3x y3z z3x z3y x2yz xy2z xyz2];
end