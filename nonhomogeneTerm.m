function symbTerm = nonhomogeneTerm(n)
    syms x y z

    terms = {};
    for order = n:-1:1 % 총 차수를 n부터 0까지 감소시키며 반복
        for max_power = order:-1:0
            for j = (order - max_power):-1:0
                k = order - max_power - j;
                terms{end+1} = x^max_power * y^j * z^k;
                terms{end+1} = y^max_power * z^j * x^k;
                terms{end+1} = z^max_power * x^j * y^k;
            end
        end
    end 
    % Convert the cell array to a symbolic array
    terms = [terms{:}];
    
    symbTerm = unique(terms,'stable');
end
