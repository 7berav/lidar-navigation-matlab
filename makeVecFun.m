function f = makeVecFun(dp, Vars)
    % dp: sym 1×nT
    % Vars: 예) {[x y z]}
    n  = numel(dp);
    fs = cell(1,n);

    for k = 1:n
        % 0 판정: simplify + isAlways 조합 (보다 강건)
        dpk_s = simplify(dp(k));
        if isequal(dpk_s, sym(0)) 
            % 영벡터를 반환하는 핸들 (N×1)
            fs{k} = @(XYZ) zeros(size(XYZ,1),1);
        else
            % 비영 항은 matlabFunction으로 일반 생성
            fs{k} = matlabFunction(dpk_s, 'Vars', Vars);
            % (필요시: 'Optimize',false 추가 가능)
        end
    end

    % 최종: N×nT 행렬로 수평 결합
    f = @(XYZ) cell2mat( cellfun(@(g) g(XYZ), fs, 'UniformOutput', false) );
end
