% Step 1: Define symbolic variables
syms x y z
order = 6;

% Step 2: Create all symbolic terms of order 6
terms = {};
for i = 0:order
    for j = 0:(order-i)
        k = order - i - j;
        terms{end+1} = x^i * y^j * z^k; % x^i * y^j * z^k where i+j+k=6
    end
end

% Convert the cell array to a symbolic array
terms = [terms{:}];

% Step 3: Define symbolic coefficients
beta = sym('beta', [1, length(terms)]); % length(terms) should match with the number of terms

% Step 4: Create the symbolic function as an inner product of beta and terms
f2 = sum(beta .* terms); % This sums up all beta(i)*term(i)

% Display the resulting symbolic function
disp(f2);