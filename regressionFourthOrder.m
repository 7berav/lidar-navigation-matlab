function output = regressionFourthOrder(inputMatrix)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    Temp=calculateFourthOrder(inputMatrix);
    
    tic
    output = inv(Temp.' * Temp)* Temp.' * ones(length(inputMatrix(:,1)),1);
    toc 


end