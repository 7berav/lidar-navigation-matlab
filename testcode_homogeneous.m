syms x y z
order = 6; % Total order for the terms

% Initialize a cell array to hold the terms
terms = {};

% Generate all combinations of exponents (i, j, k) such that i+j+k = order
for i = 0:order
    for j = 0:(order-i)
        k = order - i - j;
        terms{end+1} = x^i * y^j * z^k;
    end
end

% Convert cell array to symbolic array
symbolicArray = [terms{:}];

% Display the result
disp(symbolicArray);