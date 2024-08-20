function [beta error] = regressionFourthOrder(inputMatrix , errorInput)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    Temp=calculateDerivative(inputMatrix);
    
    tic
    beta = inv(Temp.' * Temp)* Temp.' * errorInput;
    error = ones(length(inputMatrix(:,1)),1)- Temp * beta;

    toc 
    
end