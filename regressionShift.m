function [beta error] = regressionShift(inputMatrix , errorInput)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    Temp = calculateDerivative(inputMatrix);
    
    
    beta = inv(Temp.' * Temp)* Temp.' * errorInput;
    error = errorInput- Temp * beta;

    
    
end