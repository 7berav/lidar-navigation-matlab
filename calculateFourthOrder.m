function output = calculateFourthOrder(inputMatrix,symbterm)
    % 입력 행렬에서 각 행을 추출하여 각각의 변수 x, y, z로 지정
    x1 = inputMatrix(:, 1);
    y1 = inputMatrix(:, 2);
    z1 = inputMatrix(:, 3);

    syms x y z
    numRows = size(inputMatrix, 1);
    numTerms = length(symbterm);
    
    % 출력 행렬 초기화 (n행 m열)
    output = zeros(numRows, numTerms);
    
    termFunctions = matlabFunction(symbterm, 'Vars', {[x, y, z]});
    
    % 함수 행렬을 이용하여 모든 행의 값을 한 번에 계산
    output = termFunctions(inputMatrix);
    
end
