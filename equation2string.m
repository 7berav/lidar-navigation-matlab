function poly_str = equation2string(beta , Terms)
    poly_str = "";
    if length(beta) ~= length(Terms)
        error('beta_values와 TermsB의 길이가 일치해야 합니다.');
    end
    for k = 1:length(beta)
        coeff = beta(k);
        term_k = Terms(k);
    
        % 현재 항을 문자열로 변환
        if isa(term_k, 'sym')
            term_str = char(term_k);
        elseif isstring(term_k) || ischar(term_k)
            term_str = term_k;
        else
            error('Terms의 항들은 symbolic 또는 문자열이어야 합니다.');
        end
    
        % 계수 및 항 문자열 생성
        if coeff >= 0 && k > 1
            poly_str = poly_str + " + ";
        elseif coeff < 0
            poly_str = poly_str + " - ";
            coeff = abs(coeff); % 음수 부호는 따로 처리
        end
    
        poly_str = poly_str + sprintf('%.3g*%s', coeff, term_str);
    end
    disp(poly_str);
end