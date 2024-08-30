function [beta error] = regressionShift(inputMatrix , errorInput, f,dx,dy,dz)
    %inputSensor=inputMatrix + error*rand(length(inputMatrix(:,1)),length(inputMatrix(1,:)))-error*0.5;
    
    Temp = calculateDerivative(inputMatrix, f);
    

    
    beta = inv(Temp.' * Temp)* Temp.' * errorInput;
    
    error = errorInput- Temp * beta;

    
    
end 