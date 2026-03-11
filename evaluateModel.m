function distances = evaluateModel(PP,f,beta)
    % 이 함수는 주어진 모델(beta)와 데이터(PP) 사이의 잔차를 계산합니다.
 
    distances = norm(f(PP,beta.'),1);
end