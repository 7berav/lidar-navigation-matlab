% Step 1: Define symbolic variables
syms x y z
order = 4;

% Step 2: Create all symbolic terms of order 6
terms = {};

% Loop over the maximum power of any single variable, starting from 'order' down to 0
for max_power = order:-1:0
    % For each max_power, find terms where one variable has this power,
    % and the other variables' exponents sum up to (order - max_power)
    
    % Variable 'x' has the highest power
    for j = (order-max_power):-1:0
        k = order - max_power - j;
        terms{end+1} = x^max_power * y^j * z^k;
        terms{end+1} = y^max_power * z^j * x^k;
        terms{end+1} = z^max_power * x^j * y^k;
    end
    
    %{
    % Variable 'y' has the highest power
    for k = (order-max_power):-1:0
        i = order - max_power - k;
        terms{end+1} = y^max_power * z^k * x^i;
    end
    
    % Variable 'z' has the highest power
    for i = (order-max_power):-1:0
        j = order - max_power - i;
        terms{end+1} = z^max_power * x^i * y^j;
    end
    %}
end
%%
% Convert the cell array to a symbolic array
terms = [terms{:}];
f0 = sum(terms); 
% Remove duplicate terms using unique function
uniqueTerms = unique(terms,'stable');
%uniqueTerms = unique(terms);
f05 = sum(uniqueTerms); 
% Step 3: Define symbolic coefficients
beta = sym('beta', [1, length(uniqueTerms)]); % length(terms) should match with the number of terms

%%f2 = sum(beta .* terms);

% Step 4: Create the symbolic function as an inner product of beta and terms
f1 = sum(beta .* uniqueTerms); % This sums up all beta(i)*term(i)

% Display the resulting symbolic function
disp(f1);
f1_numeric = matlabFunction(f1, 'Vars', {[x, y, z], beta});
f1_partial = @(x,y,z) f1_numeric([x,y,z], beta_values.')-1;

figure
fimplicit3(f1_partial,[-4.5 4.5 -4.5 4.5 -4.5 4.5]);
%scatter3(PPm_shift(:,1),PPm_shift(:,2),PPm_shift(:,3),'red');
axis equal
xlim([-4 4]);
ylim([-4 4]);
title('f1')