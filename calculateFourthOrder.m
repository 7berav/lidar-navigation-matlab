function output = calculateFourthOrder(inputMatrix,symbterm)
    
    if nargin(symbterm) == 1
        % 1-인자: N×3 한 번에
        output = symbterm(inputMatrix);                 % output: N×nT
    else
        % 3-인자: 열별로 나눠 전달
        output = symbterm(inputMatrix(:,1), ...
                          inputMatrix(:,2), ...
                          inputMatrix(:,3));            % output: N×nT
    end
end
