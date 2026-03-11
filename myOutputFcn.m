function stop = myOutputFcn(x, optimValues, state)
    persistent initialX initialFval finalX finalFval
    
    stop = false; % 최적화 계속 진행
    
    switch state
        case 'init'  % 최적화 시작 시점
            initialX = x;
            initialFval = optimValues.fval;

        case 'done'  % 최적화 종료 시점
            finalX = x;
            finalFval = optimValues.fval;
            
            % 최적화 종료 시마다 결과를 추가 기록
            if exist('optimization_record.mat', 'file')
                % 기존 데이터 로드 후 이어붙이기
                initialXs = [initialXs; initialX(:)'];
                initialFvals = [initialFvals; initialFval];
                finalXs = [finalXs; finalX(:)'];
                finalFvals = [finalFvals; finalFval];
            else
                % 첫번째 데이터 저장
                initialXs = initialX(:)';
                initialFvals = initialFval;
                finalXs = finalX(:)';
                finalFvals = finalFval;
            end

            % 저장
            save('optimization_record.mat', 'initialXs', 'initialFvals', 'finalXs', 'finalFvals');
    end
end