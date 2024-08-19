function [beta error] = regressionFourthOrder(inputMatrix)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    Temp=calculateFourthOrder(inputMatrix);
    
    tic
    beta = inv(Temp.' * Temp)* Temp.' * ones(length(inputMatrix(:,1)),1);
    error = ones(length(inputMatrix(:,1)),1)- Temp * beta;

    toc 
    
end