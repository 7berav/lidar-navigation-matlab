function [beta , error] = regressionShift(inputMatrix , errorInput, dxyz)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    
    %Temp = calculateDerivative(inputMatrix, f);
    

    
    beta = inv(dxyz.' * dxyz)* dxyz.' * errorInput;
    
    error = errorInput- dxyz * beta;

    
    
end 